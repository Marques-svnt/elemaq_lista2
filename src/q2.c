#include <stdio.h>
#include <math.h>
#include "q2.h" // Assumindo que q2.h declara void q2(FILE*, FILE*);

static double q2_para_radianos_static(double graus) {
    return graus * M_PI / 180.0;
}

static Q2_ForcaImperial q2_calcular_escrever_parametros_engrenagem_b(FILE *fp_respostas, double torque_lb_in) {
    Q2_ForcaImperial forca_B = {0.0, 0.0};
    fprintf(fp_respostas, "Engrenagem B (Questao 2 - recebe potencia, pinhão externo acima):\n");

    double d_gear_b_in = Q2_NUM_DENTES_GEARB / Q2_PASSO_DIAMETRAL_GEARB;
    fprintf(fp_respostas, "  Diametro Primitivo (D_gearB): %.2f in (%.4f m)\n",
            d_gear_b_in, d_gear_b_in * POLEGADA_PARA_METRO);

    double w_tb_lb = torque_lb_in / (d_gear_b_in / 2.0);
    forca_B.Fz_lb = w_tb_lb; // Horizontal (Tangencial no eixo B)
    fprintf(fp_respostas, "  Forca Tangencial no eixo devido ao pinhão (F_Bz): %.2f lb (%.2f N) (Horizontal)\n",
            forca_B.Fz_lb, forca_B.Fz_lb * LBF_PARA_NEWTON);

    double angulo_pressao_rad = q2_para_radianos_static(Q2_ANGULO_PRESSAO_GEARB_DEG);
    double w_rb_lb = w_tb_lb * tan(angulo_pressao_rad);
    forca_B.Fy_lb = -w_rb_lb; // Vertical para Baixo (Radial no eixo B)
    fprintf(fp_respostas, "  Forca Radial no eixo devido ao pinhão (F_By): %.2f lb (%.2f N) (Vertical para Baixo)\n",
            forca_B.Fy_lb, forca_B.Fy_lb * LBF_PARA_NEWTON);
    return forca_B;
}

static Q2_ForcaImperial q2_calcular_escrever_parametros_polia_d(FILE *fp_respostas, double torque_lb_in) {
    Q2_ForcaImperial forca_D = {0.0, 0.0};
    fprintf(fp_respostas, "Polia D (Questao 2 - transmite potencia):\n");
    fprintf(fp_respostas, "  Diametro da Polia D: %.2f in (%.4f m)\n",
            Q2_DIAMETRO_POLIA_D_IN, Q2_DIAMETRO_POLIA_D_IN * POLEGADA_PARA_METRO);

    double f_efetiva_polia_d_lb = torque_lb_in / (Q2_DIAMETRO_POLIA_D_IN / 2.0);
    double f_polia_d_total_lb = Q2_FATOR_FORCA_POLIA_V * f_efetiva_polia_d_lb;
    fprintf(fp_respostas, "  Forca Resultante Total na Polia (F_poliaD): %.2f lb (%.2f N) (estimada com fator %.1f)\n",
            f_polia_d_total_lb, f_polia_d_total_lb * LBF_PARA_NEWTON, Q2_FATOR_FORCA_POLIA_V);

    double angulo_forca_rad = q2_para_radianos_static(Q2_ANGULO_FORCA_POLIA_D_DEG);
    forca_D.Fy_lb = -f_polia_d_total_lb * cos(angulo_forca_rad); // Para Baixo
    forca_D.Fz_lb = f_polia_d_total_lb * sin(angulo_forca_rad);  // Para Direita
    fprintf(fp_respostas, "  Componente Vertical da Forca (F_Dy): %.2f lb (%.2f N) (Para Baixo)\n",
            forca_D.Fy_lb, forca_D.Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "  Componente Horizontal da Forca (F_Dz): %.2f lb (%.2f N) (Para Direita)\n",
            forca_D.Fz_lb, forca_D.Fz_lb * LBF_PARA_NEWTON);
    return forca_D;
}

void q2(FILE *fp_respostas, FILE *fp_dados) {
    fprintf(fp_respostas, "--- Questao 2 ---\n\n");
    fprintf(fp_respostas, "Refazendo a Questao 1 com novos parametros:\n");
    fprintf(fp_respostas, "  Potencia: %.1f HP, Rotacao: %.0f RPM\n", Q2_POTENCIA_HP, Q2_ROTACAO_RPM);
    fprintf(fp_respostas, "  Engrenagem B: %d dentes, Passo Diametral %d\n", (int)Q2_NUM_DENTES_GEARB, (int)Q2_PASSO_DIAMETRAL_GEARB);
    fprintf(fp_respostas, "  Polia D: Diametro %.1f in (%.4f m)\n\n", Q2_DIAMETRO_POLIA_D_IN, Q2_DIAMETRO_POLIA_D_IN * POLEGADA_PARA_METRO);

    // --- a) Torque (Calculado em Imperial) ---
    double torque_eixo_lb_in = (Q2_POTENCIA_HP * Q2_FATOR_TORQUE) / Q2_ROTACAO_RPM;
    fprintf(fp_respostas, "a) Magnitude do Torque no Eixo (T_BD):\n");
    fprintf(fp_respostas, "  T = (%.1f hp * %.1f / %.0f rpm) = %.2f lb.in  (%.2f N.m)\n",
            Q2_POTENCIA_HP, Q2_FATOR_TORQUE, Q2_ROTACAO_RPM, torque_eixo_lb_in, torque_eixo_lb_in * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Este torque atua entre a Engrenagem B (entrada) e a Polia D (saida).\n");
    fprintf(fp_respostas, "  Torque T_AB (antes de B) = 0 lb.in (0 N.m)\n");
    fprintf(fp_respostas, "  Torque T_DC (depois de D) = 0 lb.in (0 N.m)\n\n");

    // --- b) Forcas nos Elementos (Calculadas em Imperial, impressas com conversao) ---
    fprintf(fp_respostas, "b) Forcas nos Elementos de Transmissao:\n");
    Q2_ForcaImperial F_B = q2_calcular_escrever_parametros_engrenagem_b(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");
    Q2_ForcaImperial F_D = q2_calcular_escrever_parametros_polia_d(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");

    // --- c) Reacoes nos Mancais A e C (Calculadas em Imperial, impressas com conversao) ---
    fprintf(fp_respostas, "c) Reacoes nos Mancais (Mancal em A em x=%.1f in, Mancal em C em x=%.1f in):\n", Q2_X_A, Q2_X_C);
    double R_Ay_lb, R_Cy_lb, R_Az_lb, R_Cz_lb;

    R_Cy_lb = -(F_B.Fy_lb * Q2_X_B + F_D.Fy_lb * Q2_X_D) / Q2_X_C;
    R_Ay_lb = -(F_B.Fy_lb + F_D.Fy_lb + R_Cy_lb);
    fprintf(fp_respostas, "  Plano Vertical (XY):\n");
    fprintf(fp_respostas, "    R_Ay = %.2f lb (%.2f N)\n", R_Ay_lb, R_Ay_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    R_Cy = %.2f lb (%.2f N)\n", R_Cy_lb, R_Cy_lb * LBF_PARA_NEWTON);

    R_Cz_lb = -(F_B.Fz_lb * Q2_X_B + F_D.Fz_lb * Q2_X_D) / Q2_X_C;
    R_Az_lb = -(F_B.Fz_lb + F_D.Fz_lb + R_Cz_lb);
    fprintf(fp_respostas, "  Plano Horizontal (XZ):\n");
    fprintf(fp_respostas, "    R_Az = %.2f lb (%.2f N)\n", R_Az_lb, R_Az_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    R_Cz = %.2f lb (%.2f N)\n\n", R_Cz_lb, R_Cz_lb * LBF_PARA_NEWTON);

    // --- d) Dados para Diagramas (Convertidos para SI ao escrever em dados.txt) ---
    fprintf(fp_dados, "# Questao 2\n");
    // Formato: x[m] T[N.m] Mz[N.m] Vy[N] My[N.m] Vz[N]

    double x_in, T_lbin, Mz_lbin, Vy_lbf, My_lbin, Vz_lbf;
    double epsilon = 1e-6;

    // Ponto A (x=0)
    x_in = Q2_X_A; T_lbin = 0.0; Mz_lbin = 0.0; My_lbin = 0.0;
    Vy_lbf = 0.0; Vz_lbf = 0.0; // Antes da reacao
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);
    Vy_lbf = R_Ay_lb; Vz_lbf = R_Az_lb; // Depois da reacao
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);

    // Ponto antes de B
    x_in = Q2_X_B - epsilon; T_lbin = 0.0;
    Mz_lbin = R_Ay_lb * x_in; My_lbin = R_Az_lb * x_in;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);

    // Ponto B
    x_in = Q2_X_B;
    T_lbin = 0.0; Mz_lbin = R_Ay_lb * x_in; My_lbin = R_Az_lb * x_in; // Torque 0 antes da aplicacao
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);
    T_lbin = torque_eixo_lb_in; // Torque aplicado
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);
    Vy_lbf += F_B.Fy_lb; Vz_lbf += F_B.Fz_lb; // Forcas aplicadas
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);

    // Ponto antes de D
    x_in = Q2_X_D - epsilon; T_lbin = torque_eixo_lb_in;
    Mz_lbin = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q2_X_B);
    My_lbin = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q2_X_B);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);

    // Ponto D
    x_in = Q2_X_D;
    T_lbin = torque_eixo_lb_in; Mz_lbin = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q2_X_B); My_lbin = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q2_X_B); // Torque T_BD antes da remocao
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);
    T_lbin = 0.0; // Torque removido
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);
    Vy_lbf += F_D.Fy_lb; Vz_lbf += F_D.Fz_lb; // Forcas aplicadas
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);

    // Ponto antes de C
    x_in = Q2_X_C - epsilon; T_lbin = 0.0;
    Mz_lbin = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q2_X_B) + F_D.Fy_lb * (x_in - Q2_X_D);
    My_lbin = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q2_X_B) + F_D.Fz_lb * (x_in - Q2_X_D);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON);

    // Ponto C (Mancal)
    x_in = Q2_X_C; T_lbin = 0.0;
    Mz_lbin = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q2_X_B) + F_D.Fy_lb * (x_in - Q2_X_D);
    My_lbin = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q2_X_B) + F_D.Fz_lb * (x_in - Q2_X_D);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON); // Antes da reacao
    Vy_lbf += R_Cy_lb; Vz_lbf += R_Cz_lb; // Cortantes devem ir para ~0
    if (fabs(Vy_lbf) < epsilon) Mz_lbin = 0.0;
    if (fabs(Vz_lbf) < epsilon) My_lbin = 0.0;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin*LBIN_PARA_NM, Mz_lbin*LBIN_PARA_NM, Vy_lbf*LBF_PARA_NEWTON, My_lbin*LBIN_PARA_NM, Vz_lbf*LBF_PARA_NEWTON); // Depois da reacao

    fprintf(fp_respostas, "Dados para os diagramas (em SI) foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 2 ---\n\n\n");
}