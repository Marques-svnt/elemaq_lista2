#include <stdio.h>
#include <math.h>
#include "q3.h" // Se q3.h existir

static double q3_para_radianos_static(double graus) {
    return graus * M_PI / 180.0;
}

static double q3_calcular_torque_potencia_static(double potencia_hp, double velocidade_rpm) {
    return Q3_FATOR_TORQUE * potencia_hp / velocidade_rpm;
}

static void q3_calcular_forcas_polia_plana_A(FILE* fp_respostas, double torque_lb_in, Q3_ForcaPonto* forcaA) {
    double f1_menos_f2 = torque_lb_in / (Q3_DIAMETRO_POLIA_A_IN / 2.0);
    double f2 = f1_menos_f2 / (Q3_RELACAO_F1_F2_POLIA_A - 1.0);
    double f1 = Q3_RELACAO_F1_F2_POLIA_A * f2;
    forcaA->Fy = f1 + f2; // Para Cima
    forcaA->Fz = 0.0;
    fprintf(fp_respostas, "  Polia A (x=%.1f in, Entrada):\n", Q3_X_A);
    fprintf(fp_respostas, "    F_Ay = %.2f lb (Vertical para Cima), F_Az = %.2f lb\n", forcaA->Fy, forcaA->Fz);
}

static void q3_calcular_forcas_engrenagem_retos_C(FILE* fp_respostas, double torque_lb_in, Q3_ForcaPonto* forcaC_eng) {
    double diametro_primitivo_in = (double)Q3_NUM_DENTES_ENG_C / Q3_PASSO_DIAMETRAL_ENG_C;
    double wt = torque_lb_in / (diametro_primitivo_in / 2.0);
    double wr = wt * tan(q3_para_radianos_static(Q3_ANGULO_PRESSAO_ENG_C_DEG));
    forcaC_eng->Fy = wt; // Tangencial, para Cima no eixo Y
    forcaC_eng->Fz = wr; // Radial, Horizontal +Z
    fprintf(fp_respostas, "  Engrenagem C (x=%.1f in, Saida):\n", Q3_X_C);
    fprintf(fp_respostas, "    F_Cy (tangencial) = %.2f lb (Vertical para Cima)\n", forcaC_eng->Fy);
    fprintf(fp_respostas, "    F_Cz (radial)     = %.2f lb (Horizontal +Z)\n", forcaC_eng->Fz);
}

static void q3_calcular_forca_roda_dentada_D(FILE* fp_respostas, double torque_lb_in, Q3_ForcaPonto* forcaD_roda) {
    double fd = torque_lb_in / (Q3_DIAMETRO_RODA_DENTADA_D_IN / 2.0);
    forcaD_roda->Fy = -fd; // Para Baixo
    forcaD_roda->Fz = 0.0;
    fprintf(fp_respostas, "  Roda Dentada D (x=%.1f in, Saida):\n", Q3_X_D);
    fprintf(fp_respostas, "    F_Dy = %.2f lb (Vertical para Baixo), F_Dz = %.2f lb\n", forcaD_roda->Fy, forcaD_roda->Fz);
}

// Funcoes de calculo de diametro (apenas para respostas.txt)
static double q3_calcular_sn_linha(double sn_surf_psi, double c_size, double c_reliab){
    return sn_surf_psi * c_size * c_reliab;
}
static double q3_calcular_diametro_eixo_flexo_torcao(double M_res, double T, double Kt, double sn_linha, double sy, double N){
    if (sn_linha <= 0 || sy <= 0) return -1.0;
    double tf_sq = (Kt * M_res / sn_linha);
    tf_sq *= tf_sq;
    double tt_sq = (3.0/4.0) * (T / sy) * (T / sy);
    if (fabs(M_res) < 1e-6 && fabs(T) < 1e-6) return 0.0;
    if (fabs(M_res) < 1e-6) tf_sq = 0; // Mott p.35 (simplificado, Kt nao entra para M=0)
    if (fabs(T) < 1e-6) tt_sq = 0;
    return pow((32.0 * N / M_PI) * sqrt(tf_sq + tt_sq), 1.0/3.0);
}
static double q3_calcular_diametro_eixo_cisalhamento_puro(double V_res, double Kt, double sn_linha, double N){
    if (sn_linha <= 0) return -1.0;
    double val = (2.94 * Kt * V_res * N) / sn_linha;
    return (val < 0) ? -1.0 : sqrt(val);
}


void q3(FILE *fp_respostas, FILE *fp_dados) {
    fprintf(fp_respostas, "--- Questao 3 ---\n\n");

    // --- 1. Calculo dos Torques ---
    double t_A_entrada = q3_calcular_torque_potencia_static(Q3_POTENCIA_POLIA_A_HP, Q3_VELOCIDADE_RPM);
    double t_C_saida = q3_calcular_torque_potencia_static(Q3_POTENCIA_ENGRENAGEM_C_HP, Q3_VELOCIDADE_RPM);
    double t_D_saida = q3_calcular_torque_potencia_static(Q3_POTENCIA_RODA_DENTADA_D_HP, Q3_VELOCIDADE_RPM);

    fprintf(fp_respostas, "a) Torques nos Componentes e no Eixo:\n");
    fprintf(fp_respostas, "  Torque na Polia A (T_A_entrada): %.2f lb.in\n", t_A_entrada);
    fprintf(fp_respostas, "  Torque de Saida na Engrenagem C (T_C_saida): %.2f lb.in\n", t_C_saida);
    fprintf(fp_respostas, "  Torque de Saida na Roda Dentada D (T_D_saida): %.2f lb.in\n\n", t_D_saida);

    // Torques nos segmentos do eixo
    double t_eixo_AB = t_A_entrada; // A -> B (mancal)
    double t_eixo_BC = t_A_entrada; // B (mancal) -> C (engrenagem)
    double t_eixo_CD = t_A_entrada - t_C_saida; // C -> D (roda)
    double t_eixo_DE = t_eixo_CD - t_D_saida;   // D -> E (mancal), deve ser proximo de zero

    fprintf(fp_respostas, "  Diagrama de Torque no Eixo:\n");
    fprintf(fp_respostas, "    Segmento A (x=%.1f) a B (x=%.1f): T_AB = %.2f lb.in\n", Q3_X_A, Q3_X_B, t_eixo_AB);
    fprintf(fp_respostas, "    Segmento B (x=%.1f) a C (x=%.1f): T_BC = %.2f lb.in\n", Q3_X_B, Q3_X_C, t_eixo_BC);
    fprintf(fp_respostas, "    Segmento C (x=%.1f) a D (x=%.1f): T_CD = %.2f lb.in\n", Q3_X_C, Q3_X_D, t_eixo_CD);
    fprintf(fp_respostas, "    Segmento D (x=%.1f) a E (x=%.1f): T_DE = %.2f lb.in (deve ser ~0)\n\n", Q3_X_D, Q3_X_E, t_eixo_DE);

    // --- b) Forcas nos Elementos ---
    fprintf(fp_respostas, "b) Forcas nos Elementos de Transmissao (atuando sobre o eixo):\n");
    Q3_ForcaPonto F_A, F_C_eng, F_D_roda;
    q3_calcular_forcas_polia_plana_A(fp_respostas, t_A_entrada, &F_A);
    q3_calcular_forcas_engrenagem_retos_C(fp_respostas, t_C_saida, &F_C_eng);
    q3_calcular_forca_roda_dentada_D(fp_respostas, t_D_saida, &F_D_roda);
    fprintf(fp_respostas, "\n");

    // --- c) Reacoes nos Mancais B e E ---
    fprintf(fp_respostas, "c) Reacoes nos Mancais (Mancal B em x=%.1f, Mancal E em x=%.1f):\n", Q3_X_B, Q3_X_E);
    double R_By, R_Ey, R_Bz, R_Ez;

    // Plano Vertical (XY) - Forcas: F_A.Fy (em X_A), F_C_eng.Fy (em X_C), F_D_roda.Fy (em X_D)
    // Sum M_B_z = 0: R_Ey*(X_E-X_B) + F_A.Fy*(X_A-X_B) + F_C_eng.Fy*(X_C-X_B) + F_D_roda.Fy*(X_D-X_B) = 0
    if (fabs(Q3_X_E - Q3_X_B) < 1e-6) { // Evita divisao por zero se mancais coincidirem
        R_Ey = 0; // Ou tratar erro
    } else {
        R_Ey = -(F_A.Fy * (Q3_X_A - Q3_X_B) + F_C_eng.Fy * (Q3_X_C - Q3_X_B) + F_D_roda.Fy * (Q3_X_D - Q3_X_B)) / (Q3_X_E - Q3_X_B);
    }
    R_By = -(F_A.Fy + F_C_eng.Fy + F_D_roda.Fy + R_Ey);

    fprintf(fp_respostas, "  Plano Vertical (XY):\n");
    fprintf(fp_respostas, "    R_By = %.2f lb\n", R_By);
    fprintf(fp_respostas, "    R_Ey = %.2f lb\n", R_Ey);

    // Plano Horizontal (XZ) - Forcas: F_A.Fz (=0), F_C_eng.Fz (em X_C), F_D_roda.Fz (=0)
    // Sum M_B_y = 0: R_Ez*(X_E-X_B) + F_C_eng.Fz*(X_C-X_B) = 0
     if (fabs(Q3_X_E - Q3_X_B) < 1e-6) {
        R_Ez = 0;
    } else {
        R_Ez = -(F_C_eng.Fz * (Q3_X_C - Q3_X_B)) / (Q3_X_E - Q3_X_B);
    }
    R_Bz = -(F_C_eng.Fz + R_Ez);

    fprintf(fp_respostas, "  Plano Horizontal (XZ):\n");
    fprintf(fp_respostas, "    R_Bz = %.2f lb\n", R_Bz);
    fprintf(fp_respostas, "    R_Ez = %.2f lb\n\n", R_Ez);

    // --- d) Dados para Diagramas (escritos em dados.txt) ---
    fprintf(fp_dados, "# Questao 3\n");
    // Formato: x T Mz Vy My Vz
    double x, T_val, Mz_val, Vy_val, My_val, Vz_val;
    double epsilon = 1e-6;

    // Ponto A (x = X_A = 0.0) - Polia (Entrada de Torque e Forca F_A)
    x = Q3_X_A;
    T_val = 0; Vy_val = 0; Vz_val = 0; Mz_val = 0; My_val = 0; // Antes de A
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);
    T_val = t_eixo_AB; Vy_val = F_A.Fy; Vz_val = F_A.Fz; // Imediatamente apos A
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto antes de B (x = X_B - epsilon)
    x = Q3_X_B - epsilon;
    // T_val e Vy_val, Vz_val sao os mesmos da saida de A
    Mz_val = Vy_val * (x - Q3_X_A);
    My_val = Vz_val * (x - Q3_X_A);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto B (x = X_B) - Mancal (Reacao R_B)
    x = Q3_X_B;
    Mz_val = Vy_val * (x - Q3_X_A); // Mz em B antes da reacao
    My_val = Vz_val * (x - Q3_X_A); // My em B antes da reacao
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Antes da reacao
    Vy_val += R_By; Vz_val += R_Bz; // Apos reacao
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Apos reacao

    // Ponto antes de C (x = X_C - epsilon)
    x = Q3_X_C - epsilon;
    // T_val = t_eixo_BC (que é igual a t_eixo_AB)
    Mz_val = F_A.Fy * (x - Q3_X_A) + R_By * (x - Q3_X_B);
    My_val = F_A.Fz * (x - Q3_X_A) + R_Bz * (x - Q3_X_B);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto C (x = X_C) - Engrenagem (Saida de Torque t_C_saida, Forca F_C_eng)
    x = Q3_X_C;
    Mz_val = F_A.Fy * (x - Q3_X_A) + R_By * (x - Q3_X_B); // Mz em C antes de F_C_eng
    My_val = F_A.Fz * (x - Q3_X_A) + R_Bz * (x - Q3_X_B); // My em C antes de F_C_eng
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Antes da saida de torque/forca
    T_val = t_eixo_CD; // Torque muda
    Vy_val += F_C_eng.Fy; Vz_val += F_C_eng.Fz; // Forcas da engrenagem aplicadas
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val); // Apos saida de torque/forca

    // Ponto antes de D (x = X_D - epsilon)
    x = Q3_X_D - epsilon;
    // T_val = t_eixo_CD
    Mz_val = F_A.Fy * (x - Q3_X_A) + R_By * (x - Q3_X_B) + F_C_eng.Fy * (x - Q3_X_C);
    My_val = F_A.Fz * (x - Q3_X_A) + R_Bz * (x - Q3_X_B) + F_C_eng.Fz * (x - Q3_X_C);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto D (x = X_D) - Roda Dentada (Saida de Torque t_D_saida, Forca F_D_roda)
    x = Q3_X_D;
    Mz_val = F_A.Fy * (x - Q3_X_A) + R_By * (x - Q3_X_B) + F_C_eng.Fy * (x - Q3_X_C);
    My_val = F_A.Fz * (x - Q3_X_A) + R_Bz * (x - Q3_X_B) + F_C_eng.Fz * (x - Q3_X_C);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);
    T_val = t_eixo_DE; // Torque muda
    Vy_val += F_D_roda.Fy; Vz_val += F_D_roda.Fz; // Forcas da roda aplicadas
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto antes de E (x = X_E - epsilon)
    x = Q3_X_E - epsilon;
    // T_val = t_eixo_DE
    Mz_val = F_A.Fy * (x - Q3_X_A) + R_By * (x - Q3_X_B) + F_C_eng.Fy * (x - Q3_X_C) + F_D_roda.Fy * (x - Q3_X_D);
    My_val = F_A.Fz * (x - Q3_X_A) + R_Bz * (x - Q3_X_B) + F_C_eng.Fz * (x - Q3_X_C) + F_D_roda.Fz * (x - Q3_X_D);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);

    // Ponto E (x = X_E) - Mancal (Reacao R_E)
    x = Q3_X_E;
    Mz_val = F_A.Fy * (x - Q3_X_A) + R_By * (x - Q3_X_B) + F_C_eng.Fy * (x - Q3_X_C) + F_D_roda.Fy * (x - Q3_X_D); // Deve ser ~0
    My_val = F_A.Fz * (x - Q3_X_A) + R_Bz * (x - Q3_X_B) + F_C_eng.Fz * (x - Q3_X_C) + F_D_roda.Fz * (x - Q3_X_D); // Deve ser ~0
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);
    Vy_val += R_Ey; Vz_val += R_Ez; // Cortantes devem ir para ~0
    // Forcar Mz e My a zero se as cortantes zerarem
    if (fabs(Vy_val) < epsilon) Mz_val = 0.0;
    if (fabs(Vz_val) < epsilon) My_val = 0.0;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T_val, Mz_val, Vy_val, My_val, Vz_val);


    // --- Impressao dos calculos de diametro em respostas.txt ---
    fprintf(fp_respostas, "\nd) Calculo dos Diametros Minimos do Eixo (Resultados apenas em respostas.txt):\n");
    double c_size_default = 0.85;
    double sn_linha_default = q3_calcular_sn_linha(Q3_S_N_SURF_PSI, c_size_default, Q3_C_RELIAB_99);

    double m_res_A = sqrt(pow(0.0,2)+pow(0.0,2)); // Momento em A (inicio, antes do mancal B) e zero
    double d_A = q3_calcular_diametro_eixo_flexo_torcao(m_res_A, t_eixo_AB, Q3_KT_POLIA_A, sn_linha_default, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto A (Polia, x=%.1f): M_res=%.2f, T=%.2f, D_A = %.3f in\n", Q3_X_A, m_res_A, t_eixo_AB, d_A);

    // Recalcular Mz, My nos pontos para M_res
    double mz_val_B = F_A.Fy * (Q3_X_B - Q3_X_A);
    double my_val_B = F_A.Fz * (Q3_X_B - Q3_X_A); // F_A.Fz e 0
    double m_res_B = sqrt(pow(mz_val_B, 2) + pow(my_val_B, 2));
    double c_size_B = (d_A > 2.0 && d_A <=3.0) ? 0.82 : 0.85; // Aproximacao grosseira baseada em d_A
    double sn_linha_B = q3_calcular_sn_linha(Q3_S_N_SURF_PSI, c_size_B, Q3_C_RELIAB_99);
    double d_B = q3_calcular_diametro_eixo_flexo_torcao(m_res_B, t_eixo_BC, Q3_KT_MANCAL_B, sn_linha_B, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto B (Mancal, x=%.1f): M_res=%.2f, T=%.2f, D_B = %.3f in (C_size=%.2f)\n", Q3_X_B, m_res_B, t_eixo_BC, d_B, c_size_B);

    double mz_val_C = F_A.Fy * (Q3_X_C - Q3_X_A) + R_By * (Q3_X_C - Q3_X_B);
    double my_val_C = F_A.Fz * (Q3_X_C - Q3_X_A) + R_Bz * (Q3_X_C - Q3_X_B);
    double m_res_C = sqrt(pow(mz_val_C, 2) + pow(my_val_C, 2));
    double d_C = q3_calcular_diametro_eixo_flexo_torcao(m_res_C, t_eixo_CD, Q3_KT_ENGRENAGEM_C, sn_linha_default, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto C (Engrenagem, x=%.1f): M_res=%.2f, T=%.2f, D_C = %.3f in\n", Q3_X_C, m_res_C, t_eixo_CD, d_C);

    double mz_val_D = F_A.Fy * (Q3_X_D - Q3_X_A) + R_By * (Q3_X_D - Q3_X_B) + F_C_eng.Fy * (Q3_X_D - Q3_X_C);
    double my_val_D = F_A.Fz * (Q3_X_D - Q3_X_A) + R_Bz * (Q3_X_D - Q3_X_B) + F_C_eng.Fz * (Q3_X_D - Q3_X_C);
    double m_res_D = sqrt(pow(mz_val_D, 2) + pow(my_val_D, 2));
    double d_D = q3_calcular_diametro_eixo_flexo_torcao(m_res_D, t_eixo_DE, Q3_KT_RODA_DENTADA_D, sn_linha_default, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto D (Roda Dentada, x=%.1f): M_res=%.2f, T=%.2f, D_D = %.3f in\n", Q3_X_D, m_res_D, t_eixo_DE, d_D);
    
    double v_res_E = sqrt(pow(R_Ey,2) + pow(R_Ez,2)); // Forca cortante resultante no mancal E e a reacao
    double d_E = q3_calcular_diametro_eixo_cisalhamento_puro(fabs(v_res_E), Q3_KT_MANCAL_E, sn_linha_default, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto E (Mancal, x=%.1f): V_res=%.2f, D_E = %.3f in (por cisalhamento)\n", Q3_X_E, fabs(v_res_E), d_E);


    fprintf(fp_respostas, "\nDados para os diagramas da Questao 3 foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 3 ---\n\n\n");
}