#include <stdio.h>
#include <math.h>
#include "q4.h" // Se q4.h existir

static double q4_para_radianos_static(double graus) {
    return graus * M_PI / 180.0;
}

static double q4_calcular_torque_static(double potencia_hp, double velocidade_rpm) {
    return Q4_FATOR_TORQUE * potencia_hp / velocidade_rpm;
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
    return pow(diametro_cubo, 1.0 / 3.0);
}


void q4(FILE *fp_respostas, FILE *fp_dados) {
    fprintf(fp_respostas, "--- Questao 4 ---\n\n");

    // --- a) Magnitude do Torque no Eixo ---
    double t_saida_B = q4_calcular_torque_static(Q4_POTENCIA_ENG_B_HP, Q4_VELOCIDADE_RPM);
    double t_saida_D = q4_calcular_torque_static(Q4_POTENCIA_POLIA_D_HP, Q4_VELOCIDADE_RPM);
    double t_saida_E = q4_calcular_torque_static(Q4_POTENCIA_POLIA_E_HP, Q4_VELOCIDADE_RPM);
    double t_entrada_C = q4_calcular_torque_static(Q4_POTENCIA_CORR_C_HP, Q4_VELOCIDADE_RPM);

    fprintf(fp_respostas, "a) Magnitude do Torque no Eixo:\n");
    fprintf(fp_respostas, "  Torque Saida Engrenagem B (T_B): %.2f lb.in\n", t_saida_B);
    fprintf(fp_respostas, "  Torque Saida Polia D (T_D): %.2f lb.in\n", t_saida_D);
    fprintf(fp_respostas, "  Torque Saida Polia E (T_E): %.2f lb.in\n", t_saida_E);
    fprintf(fp_respostas, "  Torque Entrada Corrente C (T_C_entrada): %.2f lb.in\n", t_entrada_C);

    // Torques nos segmentos do eixo (C é entrada)
    double t_eixo_AB = 0.0; // Nao ha torque antes de B, que e saida alimentada por C
    double t_eixo_BC = t_saida_B; // Torque para alimentar B (vem de C)
    double t_eixo_CD = t_entrada_C - t_saida_B; // Torque remanescente em C apos alimentar B, para D e E
    double t_eixo_DE = t_eixo_CD - t_saida_D; // Torque remanescente em D apos alimentar D, para E
    double t_eixo_EF = t_eixo_DE - t_saida_E; // Deve ser ~0

    fprintf(fp_respostas, "  Torque no segmento A-B: %.2f lb.in\n", t_eixo_AB);
    fprintf(fp_respostas, "  Torque no segmento B-C: %.2f lb.in\n", t_eixo_BC);
    fprintf(fp_respostas, "  Torque no segmento C-D: %.2f lb.in\n", t_eixo_CD);
    fprintf(fp_respostas, "  Torque no segmento D-E: %.2f lb.in\n", t_eixo_DE);
    fprintf(fp_respostas, "  Torque no segmento E-F: %.2f lb.in (deve ser ~0)\n\n", t_eixo_EF);

    // --- b) Forcas nos Elementos Transmissores de Potencia ---
    fprintf(fp_respostas, "b) Forcas nos Elementos Transmissores de Potencia (atuando sobre o eixo):\n");
    Q4_ForcaPonto F_B, F_C, F_D, F_E;

    // Engrenagem B (Pinhão no eixo transmite para engrenagem ABAIXO)
    double wt_B = t_saida_B / (Q4_DIAMETRO_ENG_B_IN / 2.0);
    double wr_B = wt_B * tan(q4_para_radianos_static(Q4_ANGULO_PRESSAO_ENG_B_DEG));
    F_B.Fy = wr_B; // Radial no eixo: para CIMA (+Y)
    F_B.Fz = wt_B; // Tangencial no eixo: assumido +Z
    fprintf(fp_respostas, "  Engrenagem B (x=%.1f in):\n", Q4_X_B);
    fprintf(fp_respostas, "    Forca Vertical (F_By - Radial): %.2f lb (Para Cima)\n", F_B.Fy);
    fprintf(fp_respostas, "    Forca Horizontal (F_Bz - Tangencial): %.2f lb (Sentido +Z)\n", F_B.Fz);

    // Corrente C (Roda conjugada no lado INFERIOR ESQUERDO, 15 graus com a VERTICAL)
    // Forca da corrente no eixo C é PARA CIMA E PARA DIREITA
    double f_pull_C = t_entrada_C / (Q4_DIAMETRO_CORR_C_IN / 2.0);
    F_C.Fy = f_pull_C * cos(q4_para_radianos_static(Q4_ANGULO_FORCA_CORR_C_DEG)); // Para CIMA (+Y)
    F_C.Fz = f_pull_C * sin(q4_para_radianos_static(Q4_ANGULO_FORCA_CORR_C_DEG)); // Para DIREITA (+Z)
    fprintf(fp_respostas, "  Corrente C (x=%.1f in):\n", Q4_X_C);
    fprintf(fp_respostas, "    Forca Vertical (F_Cy): %.2f lb (Para Cima)\n", F_C.Fy);
    fprintf(fp_respostas, "    Forca Horizontal (F_Cz): %.2f lb (Para Direita, +Z)\n", F_C.Fz);

    // Polia D (Polia conjugada a 30 graus com a VERTICAL, para BAIXO e DIREITA)
    // Forca da polia no eixo D é PARA CIMA E ESQUERDA (reação) - ou enunciado se refere a força da correia.
    // O enunciado diz "transmite [...] às polias conjugadas, conforme mostrado". A figura mostra a correia puxando.
    // Então a força NO EIXO é na direção da tração da correia.
    double f_t_efetiva_D = t_saida_D / (Q4_DIAMETRO_POLIA_D_IN / 2.0);
    double f_pull_D = Q4_FATOR_FORCA_POLIA_V * f_t_efetiva_D;
    F_D.Fy = -f_pull_D * cos(q4_para_radianos_static(Q4_ANGULO_FORCA_POLIA_D_DEG)); // Para BAIXO (-Y)
    F_D.Fz = f_pull_D * sin(q4_para_radianos_static(Q4_ANGULO_FORCA_POLIA_D_DEG));  // Para DIREITA (+Z)
    fprintf(fp_respostas, "  Polia D (x=%.1f in):\n", Q4_X_D);
    fprintf(fp_respostas, "    Forca Vertical (F_Dy): %.2f lb (Para Baixo)\n", F_D.Fy);
    fprintf(fp_respostas, "    Forca Horizontal (F_Dz): %.2f lb (Para Direita, +Z)\n", F_D.Fz);

    // Polia E (Forca da polia conjugada HORIZONTAL PARA DIREITA)
    // Forca no eixo E é horizontal para DIREITA
    double f_t_efetiva_E = t_saida_E / (Q4_DIAMETRO_POLIA_E_IN / 2.0);
    double f_pull_E = Q4_FATOR_FORCA_POLIA_V * f_t_efetiva_E;
    F_E.Fy = 0.0;
    F_E.Fz = f_pull_E; // Para DIREITA (+Z)
    fprintf(fp_respostas, "  Polia E (x=%.1f in):\n", Q4_X_E);
    fprintf(fp_respostas, "    Forca Vertical (F_Ey): %.2f lb\n", F_E.Fy);
    fprintf(fp_respostas, "    Forca Horizontal (F_Ez): %.2f lb (Para Direita, +Z)\n\n", F_E.Fz);

    // --- c) Reacoes nos Rolamentos (Mancais A e F) ---
    fprintf(fp_respostas, "c) Reacoes nos Rolamentos (Mancais A em x=%.1f e F em x=%.1f):\n", Q4_X_A, Q4_X_F);
    double R_Ay, R_Fy, R_Az, R_Fz;

    // Plano Vertical (XY)
    R_Fy = -(F_B.Fy * Q4_X_B + F_C.Fy * Q4_X_C + F_D.Fy * Q4_X_D + F_E.Fy * Q4_X_E) / Q4_X_F;
    R_Ay = -(F_B.Fy + F_C.Fy + F_D.Fy + F_E.Fy + R_Fy);

    // Plano Horizontal (XZ)
    R_Fz = -(F_B.Fz * Q4_X_B + F_C.Fz * Q4_X_C + F_D.Fz * Q4_X_D + F_E.Fz * Q4_X_E) / Q4_X_F;
    R_Az = -(F_B.Fz + F_C.Fz + F_D.Fz + F_E.Fz + R_Fz);

    fprintf(fp_respostas, "  Mancal A:\n");
    fprintf(fp_respostas, "    Reacao Vertical (R_Ay): %.2f lb\n", R_Ay);
    fprintf(fp_respostas, "    Reacao Horizontal (R_Az): %.2f lb\n", R_Az);
    fprintf(fp_respostas, "  Mancal F:\n");
    fprintf(fp_respostas, "    Reacao Vertical (R_Fy): %.2f lb\n", R_Fy);
    fprintf(fp_respostas, "    Reacao Horizontal (R_Fz): %.2f lb\n\n", R_Fz);

    // --- d) Dados para Diagramas (escritos em dados.txt) ---
    fprintf(fp_dados, "# Questao 4\n");
    // Formato: x T Mz Vy My Vz
    double x, T_val, Mz_val, Vy_val, My_val, Vz_val;
    double epsilon = 1e-6;

    // Ponto A (x=0) - Mancal
    x = Q4_X_A;
    T_val = t_eixo_AB; Mz_val = 0; My_val = 0;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, 0.0, My_val, 0.0); // Antes da reacao
    Vy_val = R_Ay; Vz_val = R_Az;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Apos reacao

    // Ponto antes de B
    x = Q4_X_B - epsilon; T_val = t_eixo_AB;
    Mz_val = R_Ay * x; My_val = R_Az * x;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto B (Engrenagem)
    x = Q4_X_B; T_val = t_eixo_AB; // Torque antes de B para este ponto
    Mz_val = R_Ay * x; My_val = R_Az * x;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Antes da forca/mudanca de torque
    T_val = t_eixo_BC; // Torque muda apos B (considerando que C supre B)
    Vy_val += F_B.Fy; Vz_val += F_B.Fz;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Apos forca/mudanca de torque

    // Ponto antes de C
    x = Q4_X_C - epsilon; T_val = t_eixo_BC;
    Mz_val = R_Ay * x + F_B.Fy * (x - Q4_X_B);
    My_val = R_Az * x + F_B.Fz * (x - Q4_X_B);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto C (Corrente - Entrada de Torque)
    x = Q4_X_C; T_val = t_eixo_BC; // Torque antes da entrada de T_C_entrada
    Mz_val = R_Ay * x + F_B.Fy * (x - Q4_X_B);
    My_val = R_Az * x + F_B.Fz * (x - Q4_X_B);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Antes
    T_val = t_eixo_CD; // Torque muda apos entrada em C (e saida em B)
    Vy_val += F_C.Fy; Vz_val += F_C.Fz;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Apos

    // Ponto antes de D
    x = Q4_X_D - epsilon; T_val = t_eixo_CD;
    Mz_val = R_Ay * x + F_B.Fy * (x - Q4_X_B) + F_C.Fy * (x - Q4_X_C);
    My_val = R_Az * x + F_B.Fz * (x - Q4_X_B) + F_C.Fz * (x - Q4_X_C);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto D (Polia)
    x = Q4_X_D; T_val = t_eixo_CD;
    Mz_val = R_Ay * x + F_B.Fy * (x - Q4_X_B) + F_C.Fy * (x - Q4_X_C);
    My_val = R_Az * x + F_B.Fz * (x - Q4_X_B) + F_C.Fz * (x - Q4_X_C);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);
    T_val = t_eixo_DE;
    Vy_val += F_D.Fy; Vz_val += F_D.Fz;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto antes de E
    x = Q4_X_E - epsilon; T_val = t_eixo_DE;
    Mz_val = R_Ay * x + F_B.Fy * (x - Q4_X_B) + F_C.Fy * (x - Q4_X_C) + F_D.Fy * (x - Q4_X_D);
    My_val = R_Az * x + F_B.Fz * (x - Q4_X_B) + F_C.Fz * (x - Q4_X_C) + F_D.Fz * (x - Q4_X_D);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto E (Polia)
    x = Q4_X_E; T_val = t_eixo_DE;
    Mz_val = R_Ay * x + F_B.Fy * (x - Q4_X_B) + F_C.Fy * (x - Q4_X_C) + F_D.Fy * (x - Q4_X_D);
    My_val = R_Az * x + F_B.Fz * (x - Q4_X_B) + F_C.Fz * (x - Q4_X_C) + F_D.Fz * (x - Q4_X_D);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);
    T_val = t_eixo_EF;
    Vy_val += F_E.Fy; Vz_val += F_E.Fz;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto F (x=X_F) - Mancal
    x = Q4_X_F; T_val = t_eixo_EF; // Deveria ser zero
    Mz_val = R_Ay * x + F_B.Fy * (x - Q4_X_B) + F_C.Fy * (x - Q4_X_C) + F_D.Fy * (x - Q4_X_D) + F_E.Fy * (x - Q4_X_E); // Deveria ser ~0
    My_val = R_Az * x + F_B.Fz * (x - Q4_X_B) + F_C.Fz * (x - Q4_X_C) + F_D.Fz * (x - Q4_X_D) + F_E.Fz * (x - Q4_X_E); // Deveria ser ~0
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Antes da reacao
    Vy_val += R_Fy; Vz_val += R_Fz; // Deve ir para zero
    if (fabs(Vy_val) < epsilon) Mz_val = 0.0;
    if (fabs(Vz_val) < epsilon) My_val = 0.0;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Apos reacao

    // --- Calculo dos Diametros (para respostas.txt) ---
    fprintf(fp_respostas, "\nd) Momentos Fletores Resultantes (lb.in) e Diametros Minimos (in) para o Eixo:\n");
    double sn_linha = Q4_S_N_BASE_PSI * Q4_C_SIZE * Q4_C_RELIAB_99;
    fprintf(fp_respostas, "   (Usando Sn_linha = %.0f psi, Sy = %.0f psi, N = %.1f, Kt_anel = %.1f)\n",
           sn_linha, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N, Q4_KT_ANEL_RETENCAO);

    // Momentos e Diametros nos pontos criticos (B, C, D, E)
    double mz_calc, my_calc, mr_calc, t_calc, d_calc;

    // Ponto B
    mz_calc = R_Ay * Q4_X_B; my_calc = R_Az * Q4_X_B;
    mr_calc = sqrt(pow(mz_calc,2)+pow(my_calc,2)); t_calc = t_eixo_BC;
    d_calc = q4_calcular_diametro_eixo_ASME_static(mr_calc, t_calc, Q4_KT_ANEL_RETENCAO, sn_linha, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto B (x=%.1f): M_res = %.2f, T = %.2f, D_min = %.3f\n", Q4_X_B, mr_calc, t_calc, d_calc);

    // Ponto C
    mz_calc = R_Ay * Q4_X_C + F_B.Fy * (Q4_X_C - Q4_X_B);
    my_calc = R_Az * Q4_X_C + F_B.Fz * (Q4_X_C - Q4_X_B);
    mr_calc = sqrt(pow(mz_calc,2)+pow(my_calc,2)); t_calc = fmax(fabs(t_eixo_BC), fabs(t_eixo_CD)); // Conservador para T em C
    d_calc = q4_calcular_diametro_eixo_ASME_static(mr_calc, t_calc, Q4_KT_ANEL_RETENCAO, sn_linha, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto C (x=%.1f): M_res = %.2f, T_adj_max = %.2f, D_min = %.3f\n", Q4_X_C, mr_calc, t_calc, d_calc);

    // Ponto D
    mz_calc = R_Ay * Q4_X_D + F_B.Fy * (Q4_X_D - Q4_X_B) + F_C.Fy * (Q4_X_D - Q4_X_C);
    my_calc = R_Az * Q4_X_D + F_B.Fz * (Q4_X_D - Q4_X_B) + F_C.Fz * (Q4_X_D - Q4_X_C);
    mr_calc = sqrt(pow(mz_calc,2)+pow(my_calc,2)); t_calc = t_eixo_CD;
    d_calc = q4_calcular_diametro_eixo_ASME_static(mr_calc, t_calc, Q4_KT_ANEL_RETENCAO, sn_linha, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto D (x=%.1f): M_res = %.2f, T = %.2f, D_min = %.3f\n", Q4_X_D, mr_calc, t_calc, d_calc);

    // Ponto E
    // Calculando momentos em E a partir do mancal F (direita para esquerda)
    mz_calc = -R_Fy * (Q4_X_F - Q4_X_E); // M_E = R_Fy * (X_E - X_F)
    my_calc = -R_Fz * (Q4_X_F - Q4_X_E); // M_E = R_Fz * (X_E - X_F)
    mr_calc = sqrt(pow(mz_calc,2)+pow(my_calc,2)); t_calc = t_eixo_DE;
    d_calc = q4_calcular_diametro_eixo_ASME_static(mr_calc, t_calc, Q4_KT_ANEL_RETENCAO, sn_linha, Q4_S_Y_PSI, Q4_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto E (x=%.1f): M_res = %.2f, T = %.2f, D_min = %.3f\n", Q4_X_E, mr_calc, t_calc, d_calc);

    fprintf(fp_respostas, "\nDados para os diagramas da Questao 4 foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 4 ---\n\n\n");
}