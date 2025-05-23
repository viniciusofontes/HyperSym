function [Sv,Dv] = mSVK1(in1,in2)
%mSVK1
%    [Sv,Dv] = mSVK1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    24-Apr-2025 19:12:42

C1_1 = in1(1);
C1_2 = in1(4);
C1_3 = in1(7);
C2_2 = in1(5);
C2_3 = in1(8);
C3_3 = in1(9);
lambda = in2(:,1);
mu = in2(:,2);
t2 = C1_2.^2;
t3 = C1_3.^2;
t4 = C2_3.^2;
t5 = mu.*2.0;
t6 = C1_2.*C1_3;
t7 = C1_1.*C2_2;
t8 = C1_1.*C2_3;
t9 = C1_2.*C2_3;
t10 = C1_3.*C2_2;
t11 = C1_3.*C2_3;
t12 = C1_1.*C3_3;
t13 = C1_2.*C3_3;
t14 = C2_2.*C3_3;
t15 = C3_3.*t7;
t16 = C2_3.*t6.*2.0;
t17 = C1_1.*t4;
t18 = C2_2.*t3;
t19 = C3_3.*t2;
t20 = -t7;
t21 = -t8;
t22 = -t10;
t23 = -t12;
t24 = -t13;
t25 = -t14;
t26 = -t16;
t27 = -t15;
t28 = -t17;
t29 = -t18;
t30 = -t19;
t31 = t6+t21;
t32 = t9+t22;
t33 = t11+t24;
t34 = t2+t20;
t35 = t3+t23;
t36 = t4+t25;
t37 = t34.^2;
t38 = t35.^2;
t39 = t36.^2;
t40 = t31.^2;
t41 = t32.^2;
t42 = t33.^2;
t43 = t17+t18+t19+t26+t27;
t44 = t15+t16+t28+t29+t30;
t45 = 1.0./t43;
t47 = log(t44);
Sv = [mu.*(C1_1-1.0)+(lambda.*t36.*t45.*t47)./2.0;mu.*(C2_2-1.0)+(lambda.*t35.*t45.*t47)./2.0;mu.*(C3_3-1.0)+(lambda.*t34.*t45.*t47)./2.0;C1_2.*mu-(lambda.*t33.*t45.*t47)./2.0;C2_3.*mu-(lambda.*t31.*t45.*t47)./2.0;C1_3.*mu-(lambda.*t32.*t45.*t47)./2.0];
if nargout > 1
    t46 = t45.^2;
    t48 = t47-1.0;
    t61 = C1_1.*lambda.*t45.*t47;
    t62 = C1_2.*lambda.*t45.*t47;
    t63 = C1_3.*lambda.*t45.*t47;
    t64 = C2_2.*lambda.*t45.*t47;
    t65 = C2_3.*lambda.*t45.*t47;
    t66 = C3_3.*lambda.*t45.*t47;
    t49 = lambda.*t34.*t35.*t46;
    t50 = lambda.*t34.*t36.*t46;
    t51 = lambda.*t35.*t36.*t46;
    t52 = lambda.*t33.*t34.*t46;
    t53 = lambda.*t32.*t35.*t46;
    t54 = lambda.*t31.*t36.*t46;
    t55 = lambda.*t31.*t32.*t46;
    t56 = lambda.*t31.*t33.*t46;
    t57 = lambda.*t32.*t33.*t46;
    t67 = -t61;
    t68 = -t64;
    t69 = -t66;
    t70 = t62./2.0;
    t71 = t63./2.0;
    t72 = t65./2.0;
    t91 = lambda.*t31.*t34.*t46.*t48;
    t92 = lambda.*t32.*t34.*t46.*t48;
    t93 = lambda.*t31.*t35.*t46.*t48;
    t94 = lambda.*t33.*t35.*t46.*t48;
    t95 = lambda.*t32.*t36.*t46.*t48;
    t96 = lambda.*t33.*t36.*t46.*t48;
    t58 = -t52;
    t59 = -t53;
    t60 = -t54;
    t73 = -t70;
    t74 = -t71;
    t75 = -t72;
    t76 = t47.*t49;
    t77 = t47.*t50;
    t78 = t47.*t51;
    t79 = t47.*t52;
    t80 = t47.*t53;
    t81 = t47.*t54;
    t82 = t47.*t55;
    t83 = t47.*t56;
    t84 = t47.*t57;
    t85 = -t76;
    t86 = -t77;
    t87 = -t78;
    t88 = -t82;
    t89 = -t83;
    t90 = -t84;
    t97 = t58+t62+t79;
    t98 = t59+t63+t80;
    t99 = t60+t65+t81;
    t100 = t49+t67+t85;
    t101 = t50+t68+t86;
    t102 = t51+t69+t87;
    t103 = t55+t73+t88;
    t104 = t56+t74+t89;
    t105 = t57+t75+t90;
    Dv = [t5+lambda.*t39.*t46-lambda.*t39.*t46.*t47;t102;t101;t96;t99;t95;t102;t5+lambda.*t38.*t46-lambda.*t38.*t46.*t47;t100;t94;t93;t98;t101;t100;t5+lambda.*t37.*t46-lambda.*t37.*t46.*t47;t97;t91;t92;t96;t94;t97;mu+t66./2.0+lambda.*t42.*t46-lambda.*t42.*t46.*t47;t104;t105;t99;t93;t91;t104;mu+t61./2.0+lambda.*t40.*t46-lambda.*t40.*t46.*t47;t103;t95;t98;t92;t105;t103;mu+t64./2.0+lambda.*t41.*t46-lambda.*t41.*t46.*t47];
end
end
