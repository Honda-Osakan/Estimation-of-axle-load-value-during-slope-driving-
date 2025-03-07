# zikuzyuchisuitei
saka arxmodel

clear; clc;

% ================================
% 1. シミュレーション設定
% ================================
m1 = 1900;  m2 = 50;  m3 = 50;

k1 = 33e3;  k2 = 33e3;  % 前後サスペンション剛性
k3 = 260e3; k4 = 260e3;
c1 = 300;   c2 = 300;
l1 = 1.34;  l2 = 1.46;
J = 3000;

% 直角三角形の坂
height_slope = 30;     % 最高点 [m]
length_slope = 300;    % 水平距離 [m]

% 速度 & シミュレーション時間
v = 20;      
t_end = 10;  
dt = 0.001;  
t = 0:dt:t_end;

% 初期状態
x1 = zeros(size(t)); dx1 = zeros(size(t));
theta = zeros(size(t)); dtheta = zeros(size(t));
x2 = zeros(size(t)); dx2 = zeros(size(t));
x3 = zeros(size(t)); dx3 = zeros(size(t));
Ff = zeros(size(t)); Fr = zeros(size(t));
Fy1_array = zeros(size(t)); 
Fy2_array = zeros(size(t));

% 路面入力 (三角形坂)
road_slope = @(s) (s<0).*0 ...
                + (s>=0 & s<=length_slope).*( s/length_slope * height_slope ) ...
                + (s>length_slope).*height_slope;

% ================================
% 2. シミュレーションループ
% ================================
for i = 1:length(t)-1

    % 前後同じ路面高さとする
    road_height = road_slope(v*t(i));
    y1 = road_height;
    y2 = road_height;

    % サスペンション力
    Ff(i) = -k1 * ((x1(i) - l1*theta(i)) - x2(i)) ...
            - c1 * ((dx1(i) - l1*dtheta(i)) - dx2(i));
    Fr(i) = -k2 * ((x1(i) + l2*theta(i)) - x3(i)) ...
            - c2 * ((dx1(i) + l2*dtheta(i)) - dx3(i));

    % タイヤ力
    Fy1 = -k3*(x2(i) - y1);
    Fy2 = -k4*(x3(i) - y2);
    Fy1_array(i) = Fy1;
    Fy2_array(i) = Fy2;

    % 車体（質量 m1, 慣性 J）の運動方程式
    ddx1    = 2*(Ff(i)+Fr(i))/m1;
    ddtheta = 2*(-Ff(i)*l1 + Fr(i)*l2)/J;

    % 車軸（質量 m2, m3）の運動方程式
    ddx2 = (-k3*(x2(i)-y1) - Ff(i))/m2;
    ddx3 = (-k4*(x3(i)-y2) - Fr(i))/m3;

    % オイラー法で更新
    dx1(i+1) = dx1(i) + ddx1*dt;
    x1(i+1)  = x1(i)  + dx1(i)*dt;

    dtheta(i+1) = dtheta(i) + ddtheta*dt;
    theta(i+1)  = theta(i)  + dtheta(i)*dt;

    dx2(i+1) = dx2(i) + ddx2*dt;
    x2(i+1)  = x2(i)  + dx2(i)*dt;

    dx3(i+1) = dx3(i) + ddx3*dt;
    x3(i+1)  = x3(i)  + dx3(i)*dt;
end

% 最終ステップ補完
Ff(end) = -k1*((x1(end)-l1*theta(end)) - x2(end)) ...
          - c1*((dx1(end)-l1*dtheta(end)) - dx2(end));
Fr(end) = -k2*((x1(end)+l2*theta(end)) - x3(end)) ...
          - c2*((dx1(end)+l2*dtheta(end)) - dx3(end));

final_road_height = road_slope(v*t(end));
Fy1_array(end) = -k3*( x2(end) - final_road_height );
Fy2_array(end) = -k4*( x3(end) - final_road_height );

% ================================
% 3. タイヤ力の可視化
% ================================
figure('Color','w');
plot(t, Fy1_array, 'b','LineWidth',2); hold on;
plot(t, Fy2_array, 'r','LineWidth',2);
xlabel('時間 [s]','FontSize',16);
ylabel('垂直方向力 [N]','FontSize',16);
legend({'前車軸タイヤ力','後車軸タイヤ力'},'FontSize',16,'Location','best');
grid on; set(gca,'FontSize',16);

% ================================
% 4. ARXモデル用 データ作成 (前車軸)
% ================================
%   入力: 路面高さ (road_slope(v*t)) 
%   出力: 前車軸タイヤ力 Fy1_array (全サンプル)
u_data = arrayfun(@(tt) road_slope(v*tt), t);  % 入力列
y_data = Fy1_array;                            % 出力列

% iddata (全サンプル使用)
data_Fy1_all = iddata(y_data', u_data', dt);

% ================================
% 5. 前車軸タイヤ力のシンプルなARXモデル
% ================================
na = 5;  % 出力自己回帰次数
nb = 5;  % 入力回帰次数
nk = 1;  % 入力遅れ

model_arx_2_front = arx(data_Fy1_all, [na nb nk]);

% 比較 & Fit率取得
[y_est_all_front, fitVal_all_front, ~] = compare(data_Fy1_all, model_arx_2_front);
y_pred_all_front = y_est_all_front.OutputData;  % モデル予測

% 誤差評価
err_all_front   = y_data' - y_pred_all_front;
RMSE_all_front  = sqrt(mean(err_all_front.^2));
maxErr_all_front= max(abs(err_all_front));

disp('--- シンプルな ARXモデル (前車軸, 全データ使用) ---');
disp(model_arx_2_front);
fprintf('Fit率         = %.2f %%\n', fitVal_all_front);
fprintf('RMSE          = %.4f [N]\n', RMSE_all_front);
fprintf('最大ピーク誤差 = %.4f [N]\n', maxErr_all_front);

% ================================
% 6. 前車軸タイヤ力の出力比較プロット
% ================================
figure('Color','w');
plot(t, y_data, 'b-', 'LineWidth',1.5); hold on;
plot(t, y_pred_all_front, 'k--','LineWidth',2.5);
xlabel('時間 [s]','FontSize',16);
ylabel('垂直方向力 [N]','FontSize',16);
legend({'モデルによる実測出力','モデルによる出力予測'}, ...
       'FontSize',16,'Location','best');
grid on; set(gca,'FontSize',16);

% ================================
% 7. 後車軸タイヤ力のARXモデル用 データ作成
%    (グラフ表示はしない)
% ================================
u_data_rear = arrayfun(@(tt) road_slope(v*tt), t);  % 入力は同じ路面高さ
y_data_rear = Fy2_array;                            % 出力: 後車軸タイヤ力

data_Fy2_all = iddata(y_data_rear', u_data_rear', dt);

% 同じ次数設定で ARX モデル作成
model_arx_2_rear = arx(data_Fy2_all, [na nb nk]);

% 比較 & Fit率取得
[y_est_all_rear, fitVal_all_rear, ~] = compare(data_Fy2_all, model_arx_2_rear);
y_pred_all_rear = y_est_all_rear.OutputData;  % モデル予測

% 誤差評価
err_all_rear   = y_data_rear' - y_pred_all_rear;
RMSE_all_rear  = sqrt(mean(err_all_rear.^2));
maxErr_all_rear= max(abs(err_all_rear));

disp('--- シンプルな ARXモデル (後車軸, 全データ使用) ---');
disp(model_arx_2_rear);
fprintf('Fit率         = %.2f %%\n', fitVal_all_rear);
fprintf('RMSE          = %.4f [N]\n', RMSE_all_rear);
fprintf('最大ピーク誤差 = %.4f [N]\n', maxErr_all_rear);
