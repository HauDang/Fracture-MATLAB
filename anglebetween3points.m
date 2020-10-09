function [ angle ] = anglebetween3points( P1,P0,P2 )
P20 = [P2(1) - P0(1),P2(2) - P0(2)];
P10 = [P1(1) - P0(1),P1(2) - P0(2)];
lengthP20 = sqrt(P20(1)^2 + P20(2)^2);
lengthP10 = sqrt(P10(1)^2 + P10(2)^2);
n1 = [P20(1)/lengthP20 P20(2)/lengthP20];
n2 = [P10(1)/lengthP10 P10(2)/lengthP10];
if lengthP10*lengthP20 == 0
    angle = 0;
else
    angle = atan2(abs(n1(1)*n2(2) - n1(2)*n2(1)),n1(1)*n2(1) + n1(2)*n2(2));
end



