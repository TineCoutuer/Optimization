%% Plot bridge
figure
hold on
axis equal
grid on
disp("test")
x = [0 7 14 21 3.5 10.5 17.5];
y = [0 0 0 0 6 6 6];
% Plot nodes
plot(x, y, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8)

% Add node labels
for i = 1:4
    text(x(i), y(i)-0.7, ['P' num2str(i)], 'HorizontalAlignment', 'center');
end

% Add node labels
for i = 5:length(x)
    text(x(i), y(i)+0.7, ['P' num2str(i)], 'HorizontalAlignment', 'center');
end

% Connect P1 with P2 and P5
plot([x(1), x(2)], [y(1), y(2)], 'b-')
plot([x(1), x(5)], [y(1), y(5)], 'b-')

% Connect P2 with P5 and P6 P3
plot([x(2), x(3)], [y(2), y(3)], 'b-')
plot([x(2), x(5)], [y(2), y(5)], 'b-')
plot([x(2), x(6)], [y(2), y(6)], 'b-')

% Connect P3 with P6 and P7 and P4
plot([x(3), x(4)], [y(3), y(4)], 'b-')
plot([x(3), x(6)], [y(3), y(6)], 'b-')
plot([x(3), x(7)], [y(3), y(7)], 'b-')

plot([x(4), x(7)], [y(4), y(7)], 'b-')
plot([x(5), x(6)], [y(5), y(6)], 'b-')
plot([x(6), x(7)], [y(6), y(7)], 'b-')