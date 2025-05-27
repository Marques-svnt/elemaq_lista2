#include <stdio.h>
#include <math.h>
#include "q4.h" // Se q4.h existir

static double q4_para_radianos_static(double graus) {
    return graus * M_PI / 180.0;
}

static double q4_calcular_torque_static(double potencia_hp, double velocidade_rpm) {
    return Q4_FATOR_TORQUE * potencia_hp / velocidade_rpm; // Retorna lb.in
}

static double q4_calcular_diametro_eixo_ASME_static(double momento_fletor_lb_in, double torque_lb_in,
                                                double Kt, double sn_linha_psi, double sy_psi, double N_fator) {
    if (sn_linha_psi <= 0 || sy_psi <= 0) return -1.0;
    double termo_flexao_quad = pow((Kt * momento_fletor_lb_in / sn_linha_psi), 2.0);
    double termo_torcao_quad = (3.0 / 4.0) * pow((torque_lb_in / sy_psi), 2.0);
    if (fabs(momento_fletor_lb_in) < 1e-6 && fabs(torque_lb_in) < 1e-6) return 0.0;
    if (fabs(momento_fletor_lb_in) < 1e-6) termo_flexao_quad = 0;
    if (fabs(torque_lb_in) < 1e-6) termo_torcao_quad = 0;
    double dentro_raiz = sqrt(termo_flexao_quad + termo_torcao_quad);
    double diametro_cubo = (32.0 * N_fator / M_PI) * dentro_raiz;
    return pow(diametro_cubo, 1.0 / 3.0); // Retorna diâmetro em polegadas
}


void q4(FILE *fp_respostas, FILE *fp_dados) {
    fprintf(fp_respostas, "--- Questao 4 ---\n\n");

    // --- a) Magnitude do Torque no Eixo (Calculado em Imperial) ---
    double t_saida_B_lbin = q4_calcular_torque_static(Q4_POTENCIA_ENG_B_HP, Q4_VELOCIDADE_RPM);
    double t_saida_D_lbin = q4_calcular_torque_static(Q4_POTENCIA_POLIA_D_HP, Q4_VELOCIDADE_RPM);
    double t_saida_E_lbin = q4_calcular_torque_static(Q4_POTENCIA_POLIA_E_HP, Q4_VELOCIDADE_RPM);
    double t_entrada_C_lbin = q4_calcular_torque_static(Q4_POTENCIA_CORR_C_HP, Q4_VELOCIDADE_RPM);

    fprintf(fp_respostas, "a) Magnitude do Torque no Eixo:\n");
    fprintf(fp_respostas, "  Torque Saida Engrenagem B (T_B): %.2f lb.in (%.2f N.m)\n", t_saida_B_lbin, t_saida_B_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque Saida Polia D (T_D): %.2f lb.in (%.2f N.m)\n", t_saida_D_lbin, t_saida_D_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque Saida Polia E (T_E): %.2f lb.in (%.2f N.m)\n", t_saida_E_lbin, t_saida_E_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque Entrada Corrente C (T_C_entrada): %.2f lb.in (%.2f N.m)\n", t_entrada_C_lbin, t_entrada_C_lbin * LBIN_PARA_NM);

    double t_eixo_AB_lbin = 0.0;
    double t_eixo_BC_lbin = t_saida_B_lbin;
    double t_eixo_CD_lbin = t_entrada_C_lbin - t_saida_B_lbin;
    double t_eixo_DE_lbin = t_eixo_CD_lbin - t_saida_D_lbin;
    double t_eixo_EF_lbin = t_eixo_DE_lbin - t_saida_E_lbin;

    fprintf(fp_respostas, "  Torque no segmento A-B (x=%.1f a %.1f in): %.2f lb.in (%.2f N.m)\n", Q4_X_A, Q4_X_B, t_eixo_AB_lbin, t_eixo_AB_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque no segmento B-C (x=%.1f a %.1f in): %.2f lb.in (%.2f N.m)\n", Q4_X_B, Q4_X_C, t_eixo_BC_lbin, t_eixo_BC_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque no segmento C-D (x=%.1f a %.1f in): %.2f lb.in (%.2f N.m)\n", Q4_X_C, Q4_X_D, t_eixo_CD_lbin, t_eixo_CD_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque no segmento D-E (x=%.1f a %.1f in): %.2f lb.in (%.2f N.m)\n", Q4_X_D, Q4_X_E, t_eixo_DE_lbin, t_eixo_DE_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque no segmento E-F (x=%.1f a %.1f in): %.2f lb.in (%.2f N.m) (deve ser ~0)\n\n", Q4_X_E, Q4_X_F, t_eixo_EF_lbin, t_eixo_EF_lbin * LBIN_PARA_NM);

    // --- b) Forcas nos Elementos Transmissores de Potencia (Imperial, impressas com conversao) ---
    fprintf(fp_respostas, "b) Forcas nos Elementos Transmissores de Potencia (atuando sobre o eixo):\n");
    Q4_ForcaImperial F_B, F_C, F_D, F_E;

    double wt_B_lb = t_saida_B_lbin / (Q4_DIAMETRO_ENG_B_IN / 2.0);
    double wr_B_lb = wt_B_lb * tan(q4_para_radianos_static(Q4_ANGULO_PRESSAO_ENG_B_DEG));
    F_B.Fy_lb = wr_B_lb; F_B.Fz_lb = wt_B_lb;
    fprintf(fp_respostas, "  Engrenagem B (x=%.1f in / %.3f m):\n", Q4_X_B, Q4_X_B * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    Forca Vertical (F_By - Radial): %.2f lb (%.2f N) (Para Cima)\n", F_B.Fy_lb, F_B.Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    Forca Horizontal (F_Bz - Tangencial): %.2f lb (%.2f N) (Sentido +Z)\n", F_B.Fz_lb, F_B.Fz_lb * LBF_PARA_NEWTON);

    double f_pull_C_lb = t_entrada_C_lbin / (Q4_DIAMETRO_CORR_C_IN / 2.0);
    F_C.Fy_lb = -f_pull_C_lb * cos(q4_para_radianos_static(Q4_ANGULO_FORCA_CORR_C_DEG));
    F_C.Fz_lb = -f_pull_C_lb * sin(q4_para_radianos_static(Q4_ANGULO_FORCA_CORR_C_DEG));
    fprintf(fp_respostas, "  Corrente C (x=%.1f in / %.3f m):\n", Q4_X_C, Q4_X_C * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    Forca Vertical (F_Cy): %.2f lb (%.2f N) (Para Baixo)\n", F_C.Fy_lb, F_C.Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    Forca Horizontal (F_Cz): %.2f lb (%.2f N) (Para Esquerda, -Z)\n", F_C.Fz_lb, F_C.Fz_lb * LBF_PARA_NEWTON);

    double f_t_efetiva_D_lb = t_saida_D_lbin / (Q4_DIAMETRO_POLIA_D_IN / 2.0);
    double f_pull_D_lb = Q4_FATOR_FORCA_POLIA_V * f_t_efetiva_D_lb;
    F_D.Fy_lb = -f_pull_D_lb * cos(q4_para_radianos_static(Q4_ANGULO_FORCA_POLIA_D_DEG));
    F_D.Fz_lb = f_pull_D_lb * sin(q4_para_radianos_static(Q4_ANGULO_FORCA_POLIA_D_DEG));
    fprintf(fp_respostas, "  Polia D (x=%.1f in / %.3f m):\n", Q4_X_D, Q4_X_D * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    Forca Vertical (F_Dy): %.2f lb (%.2f N) (Para Baixo)\n", F_D.Fy_lb, F_D.Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    Forca Horizontal (F_Dz): %.2f lb (%.2f N) (Para Direita, +Z)\n", F_D.Fz_lb, F_D.Fz_lb * LBF_PARA_NEWTON);

    double f_t_efetiva_E_lb = t_saida_E_lbin / (Q4_DIAMETRO_POLIA_E_IN / 2.0);
    double f_pull_E_lb = Q4_FATOR_FORCA_POLIA_V * f_t_efetiva_E_lb;
    F_E.Fy_lb = 0.0; F_E.Fz_lb = f_pull_E_lb;
    fprintf(fp_respostas, "  Polia E (x=%.1f in / %.3f m):\n", Q4_X_E, Q4_X_E * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    Forca Vertical (F_Ey): %.2f lb (%.2f N)\n", F_E.Fy_lb, F_E.Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    Forca Horizontal (F_Ez): %.2f lb (%.2f N) (Para Direita, +Z)\n\n", F_E.Fz_lb, F_E.Fz_lb * LBF_PARA_NEWTON);

    // --- c) Reacoes nos Rolamentos (Mancais A e F) (Imperial, impressas com conversao) ---
    fprintf(fp_respostas, "c) Reacoes nos Rolamentos (Mancais A em x=%.1f in, F em x=%.1f in):\n", Q4_X_A, Q4_X_F);
    double R_Ay_lb, R_Fy_lb, R_Az_lb, R_Fz_lb;

    R_Fy_lb = -(F_B.Fy_lb * Q4_X_B + F_C.Fy_lb * Q4_X_C + F_D.Fy_lb * Q4_X_D + F_E.Fy_lb * Q4_X_E) / Q4_X_F;
    R_Ay_lb = -(F_B.Fy_lb + F_C.Fy_lb + F_D.Fy_lb + F_E.Fy_lb + R_Fy_lb);
    R_Fz_lb = -(F_B.Fz_lb * Q4_X_B + F_C.Fz_lb * Q4_X_C + F_D.Fz_lb * Q4_X_D + F_E.Fz_lb * Q4_X_E) / Q4_X_F;
    R_Az_lb = -(F_B.Fz_lb + F_C.Fz_lb + F_D.Fz_lb + F_E.Fz_lb + R_Fz_lb);

    fprintf(fp_respostas, "  Mancal A:\n");
    fprintf(fp_respostas, "    Reacao Vertical (R_Ay): %.2f lb (%.2f N)\n", R_Ay_lb, R_Ay_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    Reacao Horizontal (R_Az): %.2f lb (%.2f N)\n", R_Az_lb, R_Az_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "  Mancal F:\n");
    fprintf(fp_respostas, "    Reacao Vertical (R_Fy): %.2f lb (%.2f N)\n", R_Fy_lb, R_Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    Reacao Horizontal (R_Fz): %.2f lb (%.2f N)\n\n", R_Fz_lb, R_Fz_lb * LBF_PARA_NEWTON);

    // --- d) Dados para Diagramas (Convertidos para SI ao escrever em dados.txt) ---
    fprintf(fp_dados, "# Questao 4\n");
    // Formato: x[m] T[N.m] Mz[N.m] Vy[N] My[N.m] Vz[N]
    double x_in, T_lbin_val, Mz_lbin_val, Vy_lbf_val, My_lbin_val, Vz_lbf_val; // <<--- CORREÇÃO AQUI
    double epsilon = 1e-6;

    // Ponto A (x=0) - Mancal
    x_in = Q4_X_A; T_lbin_val = t_eixo_AB_lbin; Mz_lbin_val = 0.0; My_lbin_val = 0.0;
    Vy_lbf_val = 0.0; Vz_lbf_val = 0.0;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    Vy_lbf_val = R_Ay_lb; Vz_lbf_val = R_Az_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto antes de B
    x_in = Q4_X_B - epsilon; T_lbin_val = t_eixo_AB_lbin;
    Mz_lbin_val = R_Ay_lb * x_in; My_lbin_val = R_Az_lb * x_in;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto B (Engrenagem)
    x_in = Q4_X_B; T_lbin_val = t_eixo_AB_lbin;
    Mz_lbin_val = R_Ay_lb * x_in; My_lbin_val = R_Az_lb * x_in;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    T_lbin_val = t_eixo_BC_lbin; Vy_lbf_val += F_B.Fy_lb; Vz_lbf_val += F_B.Fz_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto antes de C
    x_in = Q4_X_C - epsilon; T_lbin_val = t_eixo_BC_lbin;
    Mz_lbin_val = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q4_X_B);
    My_lbin_val = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q4_X_B);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto C (Corrente)
    x_in = Q4_X_C; T_lbin_val = t_eixo_BC_lbin;
    Mz_lbin_val = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q4_X_B);
    My_lbin_val = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q4_X_B);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    T_lbin_val = t_eixo_CD_lbin; Vy_lbf_val += F_C.Fy_lb; Vz_lbf_val += F_C.Fz_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto antes de D
    x_in = Q4_X_D - epsilon; T_lbin_val = t_eixo_CD_lbin;
    Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q4_X_B) + F_C.Fy_lb*(x_in-Q4_X_C);
    My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q4_X_B) + F_C.Fz_lb*(x_in-Q4_X_C);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto D (Polia)
    x_in = Q4_X_D; T_lbin_val = t_eixo_CD_lbin;
    Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q4_X_B) + F_C.Fy_lb*(x_in-Q4_X_C);
    My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q4_X_B) + F_C.Fz_lb*(x_in-Q4_X_C);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    T_lbin_val = t_eixo_DE_lbin; Vy_lbf_val += F_D.Fy_lb; Vz_lbf_val += F_D.Fz_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto antes de E
    x_in = Q4_X_E - epsilon; T_lbin_val = t_eixo_DE_lbin;
    Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q4_X_B) + F_C.Fy_lb*(x_in-Q4_X_C) + F_D.Fy_lb*(x_in-Q4_X_D);
    My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q4_X_B) + F_C.Fz_lb*(x_in-Q4_X_C) + F_D.Fz_lb*(x_in-Q4_X_D);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto E (Polia)
    x_in = Q4_X_E; T_lbin_val = t_eixo_DE_lbin;
    Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q4_X_B) + F_C.Fy_lb*(x_in-Q4_X_C) + F_D.Fy_lb*(x_in-Q4_X_D);
    My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q4_X_B) + F_C.Fz_lb*(x_in-Q4_X_C) + F_D.Fz_lb*(x_in-Q4_X_D);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    T_lbin_val = t_eixo_EF_lbin; Vy_lbf_val += F_E.Fy_lb; Vz_lbf_val += F_E.Fz_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto F (Mancal)
    x_in = Q4_X_F; T_lbin_val = t_eixo_EF_lbin; // Deve ser zero
    Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q4_X_B) + F_C.Fy_lb*(x_in-Q4_X_C) + F_D.Fy_lb*(x_in-Q4_X_D) + F_E.Fy_lb*(x_in-Q4_X_E); // ~0
    My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q4_X_B) + F_C.Fz_lb*(x_in-Q4_X_C) + F_D.Fz_lb*(x_in-Q4_X_D) + F_E.Fz_lb*(x_in-Q4_X_E); // ~0
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    Vy_lbf_val += R_Fy_lb; Vz_lbf_val += R_Fz_lb; // Devem zerar
    if (fabs(Vy_lbf_val) < epsilon) Mz_lbin_val = 0.0;
    if (fabs(Vz_lbf_val) < epsilon) My_lbin_val = 0.0;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // --- Calculo dos Diametros (Imperial, impressas com conversao para respostas.txt) ---
    fprintf(fp_respostas, "\nd) Momentos Fletores Resultantes e Diametros Minimos para o Eixo:\n");
    double sn_linha_psi = Q4_S_N_BASE_PSI * Q4_C_SIZE * Q4_C_RELIAB_99;
    fprintf(fp_respostas, "   (Usando Sn_linha = %.0f psi (%.2f MPa), Sy = %.0f psi (%.2f MPa), N = %.1f, Kt_anel = %.1f)\n",
           sn_linha_psi, sn_linha_psi * PSI_PARA_PASCAL * PASCAL_PARA_MPA,
           Q4_S_Y_PSI, Q4_S_Y_PSI * PSI_PARA_PASCAL * PASCAL_PARA_MPA,
           Q4_FATOR_PROJETO_N, Q4_KT_ANEL_RETENCAO);

    double mz_calc_lbin, my_calc_lbin, mr_calc_lbin, t_calc_lbin, d_calc_in;

    // Ponto B
    mz_calc_lbin = R_Ay_lb * Q4_X_B; my_calc_lbin = R_Az_lb * Q4_X_B;
    mr_calc_lbin = sqrt(pow(mz_calc_lbin,2)+pow(my_calc_lbin,2)); t_calc_lbin = t_eixo_BC_lbin;
    d_calc_in = q4_calcular_diametro_eixo_ASME_static(mr_calc_lbin, t_calc_lbin, Q4_KT_ANEL_RETENCAO, sn_linha_psi, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto B (x=%.1f in): M_res = %.2f lb.in (%.2f N.m), T = %.2f lb.in (%.2f N.m), D_min = %.3f in (%.2f mm)\n", Q4_X_B, mr_calc_lbin, mr_calc_lbin*LBIN_PARA_NM, t_calc_lbin, t_calc_lbin*LBIN_PARA_NM, d_calc_in, d_calc_in*POLEGADA_PARA_METRO*1000);

    // Ponto C
    mz_calc_lbin = R_Ay_lb * Q4_X_C + F_B.Fy_lb * (Q4_X_C - Q4_X_B);
    my_calc_lbin = R_Az_lb * Q4_X_C + F_B.Fz_lb * (Q4_X_C - Q4_X_B);
    mr_calc_lbin = sqrt(pow(mz_calc_lbin,2)+pow(my_calc_lbin,2)); t_calc_lbin = fmax(fabs(t_eixo_BC_lbin), fabs(t_eixo_CD_lbin));
    d_calc_in = q4_calcular_diametro_eixo_ASME_static(mr_calc_lbin, t_calc_lbin, Q4_KT_ANEL_RETENCAO, sn_linha_psi, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto C (x=%.1f in): M_res = %.2f lb.in (%.2f N.m), T_adj_max = %.2f lb.in (%.2f N.m), D_min = %.3f in (%.2f mm)\n", Q4_X_C, mr_calc_lbin, mr_calc_lbin*LBIN_PARA_NM, t_calc_lbin, t_calc_lbin*LBIN_PARA_NM, d_calc_in, d_calc_in*POLEGADA_PARA_METRO*1000);

    // Ponto D
    mz_calc_lbin = R_Ay_lb*Q4_X_D + F_B.Fy_lb*(Q4_X_D-Q4_X_B) + F_C.Fy_lb*(Q4_X_D-Q4_X_C);
    my_calc_lbin = R_Az_lb*Q4_X_D + F_B.Fz_lb*(Q4_X_D-Q4_X_B) + F_C.Fz_lb*(Q4_X_D-Q4_X_C);
    mr_calc_lbin = sqrt(pow(mz_calc_lbin,2)+pow(my_calc_lbin,2)); t_calc_lbin = t_eixo_DE_lbin; // Torque no segmento DE (após D)
    d_calc_in = q4_calcular_diametro_eixo_ASME_static(mr_calc_lbin, t_calc_lbin, Q4_KT_ANEL_RETENCAO, sn_linha_psi, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto D (x=%.1f in): M_res = %.2f lb.in (%.2f N.m), T = %.2f lb.in (%.2f N.m), D_min = %.3f in (%.2f mm)\n", Q4_X_D, mr_calc_lbin, mr_calc_lbin*LBIN_PARA_NM, t_calc_lbin, t_calc_lbin*LBIN_PARA_NM, d_calc_in, d_calc_in*POLEGADA_PARA_METRO*1000);

    // Ponto E
    mz_calc_lbin = -R_Fy_lb * (Q4_X_F - Q4_X_E);
    my_calc_lbin = -R_Fz_lb * (Q4_X_F - Q4_X_E);
    mr_calc_lbin = sqrt(pow(mz_calc_lbin,2)+pow(my_calc_lbin,2)); t_calc_lbin = t_eixo_EF_lbin; // Torque no segmento EF (após E, deve ser 0)
    d_calc_in = q4_calcular_diametro_eixo_ASME_static(mr_calc_lbin, t_calc_lbin, Q4_KT_ANEL_RETENCAO, sn_linha_psi, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto E (x=%.1f in): M_res = %.2f lb.in (%.2f N.m), T = %.2f lb.in (%.2f N.m), D_min = %.3f in (%.2f mm)\n", Q4_X_E, mr_calc_lbin, mr_calc_lbin*LBIN_PARA_NM, t_calc_lbin, t_calc_lbin*LBIN_PARA_NM, d_calc_in, d_calc_in*POLEGADA_PARA_METRO*1000);

    fprintf(fp_respostas, "\nDados para os diagramas (em SI) foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 4 ---\n\n\n");
}