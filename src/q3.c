#include <stdio.h>
#include <math.h>
#include "q3.h" // Se q3.h existir

static double q3_para_radianos_static(double graus)
{
    return graus * M_PI / 180.0;
}

static double q3_calcular_torque_potencia_static(double potencia_hp, double velocidade_rpm)
{
    return Q3_FATOR_TORQUE * potencia_hp / velocidade_rpm;
}

static void q3_calcular_forcas_polia_plana_A(FILE *fp_respostas, double torque_lb_in, Q3_ForcaImperial *forcaA)
{
    double f1_menos_f2 = torque_lb_in / (Q3_DIAMETRO_POLIA_A_IN / 2.0);
    double f2 = f1_menos_f2 / (Q3_RELACAO_F1_F2_POLIA_A - 1.0);
    double f1 = Q3_RELACAO_F1_F2_POLIA_A * f2;
    forcaA->Fy_lb = f1 + f2; // Para Cima
    forcaA->Fz_lb = 0.0;

    fprintf(fp_respostas, "  Polia A (x=%.1f in / %.3f m, Entrada):\n", Q3_X_A, Q3_X_A * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    Diametro: %.2f in (%.4f m)\n", Q3_DIAMETRO_POLIA_A_IN, Q3_DIAMETRO_POLIA_A_IN * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    F_Ay = %.2f lb (%.2f N) (Vertical para Cima)\n", forcaA->Fy_lb, forcaA->Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    F_Az = %.2f lb (%.2f N)\n", forcaA->Fz_lb, forcaA->Fz_lb * LBF_PARA_NEWTON);
}

static void q3_calcular_forcas_engrenagem_retos_C(FILE *fp_respostas, double torque_lb_in, Q3_ForcaImperial *forcaC_eng)
{
    double diametro_primitivo_in = (double)Q3_NUM_DENTES_ENG_C / Q3_PASSO_DIAMETRAL_ENG_C;
    double wt_lb = torque_lb_in / (diametro_primitivo_in / 2.0);
    double wr_lb = wt_lb * tan(q3_para_radianos_static(Q3_ANGULO_PRESSAO_ENG_C_DEG));

    forcaC_eng->Fy_lb = wt_lb; // Tangencial, para Cima no eixo Y (conforme código original)
    forcaC_eng->Fz_lb = wr_lb; // Radial, Horizontal +Z (conforme código original)

    fprintf(fp_respostas, "  Engrenagem C (x=%.1f in / %.3f m, Saida):\n", Q3_X_C, Q3_X_C * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    Diametro Primitivo: %.2f in (%.4f m)\n", diametro_primitivo_in, diametro_primitivo_in * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    F_Cy (tangencial) = %.2f lb (%.2f N) (Vertical para Cima)\n", forcaC_eng->Fy_lb, forcaC_eng->Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    F_Cz (radial)     = %.2f lb (%.2f N) (Horizontal +Z)\n", forcaC_eng->Fz_lb, forcaC_eng->Fz_lb * LBF_PARA_NEWTON);
}

static void q3_calcular_forca_roda_dentada_D(FILE *fp_respostas, double torque_lb_in, Q3_ForcaImperial *forcaD_roda)
{
    double fd_lb = torque_lb_in / (Q3_DIAMETRO_RODA_DENTADA_D_IN / 2.0);
    forcaD_roda->Fy_lb = -fd_lb; // Para Baixo
    forcaD_roda->Fz_lb = 0.0;

    fprintf(fp_respostas, "  Roda Dentada D (x=%.1f in / %.3f m, Saida):\n", Q3_X_D, Q3_X_D * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    Diametro: %.2f in (%.4f m)\n", Q3_DIAMETRO_RODA_DENTADA_D_IN, Q3_DIAMETRO_RODA_DENTADA_D_IN * POLEGADA_PARA_METRO);
    fprintf(fp_respostas, "    F_Dy = %.2f lb (%.2f N) (Vertical para Baixo)\n", forcaD_roda->Fy_lb, forcaD_roda->Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    F_Dz = %.2f lb (%.2f N)\n", forcaD_roda->Fz_lb, forcaD_roda->Fz_lb * LBF_PARA_NEWTON);
}

static double q3_calcular_sn_linha(double sn_surf_psi, double c_size, double c_reliab)
{
    return sn_surf_psi * c_size * c_reliab;
}
static double q3_calcular_diametro_eixo_flexo_torcao(double M_res_lbin, double T_lbin, double Kt,
                                                     double sn_linha_psi, double sy_psi, double N_fator)
{
    if (sn_linha_psi <= 0 || sy_psi <= 0)
        return -1.0;
    double tf_sq = pow((Kt * M_res_lbin / sn_linha_psi), 2.0);
    double tt_sq = (3.0 / 4.0) * pow((T_lbin / sy_psi), 2.0);
    if (fabs(M_res_lbin) < 1e-6 && fabs(T_lbin) < 1e-6)
        return 0.0;
    if (fabs(M_res_lbin) < 1e-6)
        tf_sq = 0;
    if (fabs(T_lbin) < 1e-6)
        tt_sq = 0;
    double dentro_raiz = sqrt(tf_sq + tt_sq);
    double diametro_cubo = (32.0 * N_fator / M_PI) * dentro_raiz;
    return pow(diametro_cubo, 1.0 / 3.0); // Retorna diâmetro em polegadas
}
static double q3_calcular_diametro_eixo_cisalhamento_puro(double V_res_lbf, double Kt,
                                                          double sn_linha_psi, double N_fator)
{
    if (sn_linha_psi <= 0)
        return -1.0;
    double val = (2.94 * Kt * V_res_lbf * N_fator) / sn_linha_psi;
    return (val < 0) ? -1.0 : sqrt(val); // Retorna diâmetro em polegadas
}

void q3(FILE *fp_respostas, FILE *fp_dados)
{
    fprintf(fp_respostas, "--- Questao 3 ---\n\n");

    // --- 1. Calculo dos Torques (Imperial) ---
    double t_A_entrada_lbin = q3_calcular_torque_potencia_static(Q3_POTENCIA_POLIA_A_HP, Q3_VELOCIDADE_RPM);
    double t_C_saida_lbin = q3_calcular_torque_potencia_static(Q3_POTENCIA_ENGRENAGEM_C_HP, Q3_VELOCIDADE_RPM);
    double t_D_saida_lbin = q3_calcular_torque_potencia_static(Q3_POTENCIA_RODA_DENTADA_D_HP, Q3_VELOCIDADE_RPM);

    fprintf(fp_respostas, "a) Torques nos Componentes e no Eixo:\n");
    fprintf(fp_respostas, "  Torque na Polia A (T_A_entrada): %.2f lb.in (%.2f N.m)\n", t_A_entrada_lbin, t_A_entrada_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque de Saida na Engrenagem C (T_C_saida): %.2f lb.in (%.2f N.m)\n", t_C_saida_lbin, t_C_saida_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Torque de Saida na Roda Dentada D (T_D_saida): %.2f lb.in (%.2f N.m)\n\n", t_D_saida_lbin, t_D_saida_lbin * LBIN_PARA_NM);

    double t_eixo_AB_lbin = t_A_entrada_lbin;
    double t_eixo_BC_lbin = t_A_entrada_lbin;
    double t_eixo_CD_lbin = t_A_entrada_lbin - t_C_saida_lbin;
    double t_eixo_DE_lbin = t_eixo_CD_lbin - t_D_saida_lbin;

    fprintf(fp_respostas, "  Diagrama de Torque no Eixo:\n");
    fprintf(fp_respostas, "    Segmento A (x=%.1f) a B (x=%.1f): T_AB = %.2f lb.in (%.2f N.m)\n", Q3_X_A, Q3_X_B, t_eixo_AB_lbin, t_eixo_AB_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "    Segmento B (x=%.1f) a C (x=%.1f): T_BC = %.2f lb.in (%.2f N.m)\n", Q3_X_B, Q3_X_C, t_eixo_BC_lbin, t_eixo_BC_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "    Segmento C (x=%.1f) a D (x=%.1f): T_CD = %.2f lb.in (%.2f N.m)\n", Q3_X_C, Q3_X_D, t_eixo_CD_lbin, t_eixo_CD_lbin * LBIN_PARA_NM);
    fprintf(fp_respostas, "    Segmento D (x=%.1f) a E (x=%.1f): T_DE = %.2f lb.in (%.2f N.m) (deve ser ~0)\n\n", Q3_X_D, Q3_X_E, t_eixo_DE_lbin, t_eixo_DE_lbin * LBIN_PARA_NM);

    // --- b) Forcas nos Elementos (Imperial, impressas com conversao) ---
    fprintf(fp_respostas, "b) Forcas nos Elementos de Transmissao (atuando sobre o eixo):\n");
    Q3_ForcaImperial F_A, F_C_eng, F_D_roda;
    q3_calcular_forcas_polia_plana_A(fp_respostas, t_A_entrada_lbin, &F_A);
    q3_calcular_forcas_engrenagem_retos_C(fp_respostas, t_C_saida_lbin, &F_C_eng);
    q3_calcular_forca_roda_dentada_D(fp_respostas, t_D_saida_lbin, &F_D_roda);
    fprintf(fp_respostas, "\n");

    // --- c) Reacoes nos Mancais B e E (Imperial, impressas com conversao) ---
    fprintf(fp_respostas, "c) Reacoes nos Mancais (Mancal B em x=%.1f in, Mancal E em x=%.1f in):\n", Q3_X_B, Q3_X_E);
    double R_By_lb, R_Ey_lb, R_Bz_lb, R_Ez_lb;

    if (fabs(Q3_X_E - Q3_X_B) < 1e-6)
    {
        R_Ey_lb = 0;
        R_Ez_lb = 0;
    } // Evita divisao por zero
    else
    {
        R_Ey_lb = -(F_A.Fy_lb * (Q3_X_A - Q3_X_B) + F_C_eng.Fy_lb * (Q3_X_C - Q3_X_B) + F_D_roda.Fy_lb * (Q3_X_D - Q3_X_B)) / (Q3_X_E - Q3_X_B);
        R_Ez_lb = -(F_A.Fz_lb * (Q3_X_A - Q3_X_B) + F_C_eng.Fz_lb * (Q3_X_C - Q3_X_B) + F_D_roda.Fz_lb * (Q3_X_D - Q3_X_B)) / (Q3_X_E - Q3_X_B);
        // Nota: F_A.Fz_lb e F_D_roda.Fz_lb sao 0 neste problema.
    }
    R_By_lb = -(F_A.Fy_lb + F_C_eng.Fy_lb + F_D_roda.Fy_lb + R_Ey_lb);
    R_Bz_lb = -(F_A.Fz_lb + F_C_eng.Fz_lb + F_D_roda.Fz_lb + R_Ez_lb);

    fprintf(fp_respostas, "  Plano Vertical (XY):\n");
    fprintf(fp_respostas, "    R_By = %.2f lb (%.2f N)\n", R_By_lb, R_By_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    R_Ey = %.2f lb (%.2f N)\n", R_Ey_lb, R_Ey_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "  Plano Horizontal (XZ):\n");
    fprintf(fp_respostas, "    R_Bz = %.2f lb (%.2f N)\n", R_Bz_lb, R_Bz_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    R_Ez = %.2f lb (%.2f N)\n\n", R_Ez_lb, R_Ez_lb * LBF_PARA_NEWTON);

    // --- d) Dados para Diagramas (Convertidos para SI ao escrever em dados.txt) ---
    fprintf(fp_dados, "# Questao 3\n");
    // Formato: x[m] T[N.m] Mz[N.m] Vy[N] My[N.m] Vz[N]
    double x_in, T_lbin_val, Mz_lbin_val, Vy_lbf_val, My_lbin_val, Vz_lbf_val;
    double epsilon = 1e-6;

    // Ponto A (x = X_A)
    x_in = Q3_X_A;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, 0.0 * LBIN_PARA_NM, 0.0 * LBIN_PARA_NM, 0.0 * LBF_PARA_NEWTON, 0.0 * LBIN_PARA_NM, 0.0 * LBF_PARA_NEWTON); // Antes de F_A
    T_lbin_val = t_eixo_AB_lbin;
    Vy_lbf_val = F_A.Fy_lb;
    Vz_lbf_val = F_A.Fz_lb;
    Mz_lbin_val = 0.0;
    My_lbin_val = 0.0;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

    // Ponto antes de B
    x_in = Q3_X_B - epsilon;
    T_lbin_val = t_eixo_AB_lbin;
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A);
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

    // Ponto B (Mancal)
    x_in = Q3_X_B;
    T_lbin_val = t_eixo_BC_lbin; // Torque continua t_eixo_AB = t_eixo_BC
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A);
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Antes da reacao
    Vy_lbf_val += R_By_lb;
    Vz_lbf_val += R_Bz_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Apos reacao

    // Ponto antes de C
    x_in = Q3_X_C - epsilon;
    T_lbin_val = t_eixo_BC_lbin;
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A) + R_By_lb * (x_in - Q3_X_B);
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A) + R_Bz_lb * (x_in - Q3_X_B);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

    // Ponto C (Engrenagem)
    x_in = Q3_X_C;
    T_lbin_val = t_eixo_BC_lbin;
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A) + R_By_lb * (x_in - Q3_X_B);
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A) + R_Bz_lb * (x_in - Q3_X_B);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Antes
    T_lbin_val = t_eixo_CD_lbin;
    Vy_lbf_val += F_C_eng.Fy_lb;
    Vz_lbf_val += F_C_eng.Fz_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Apos

    // Ponto antes de D
    x_in = Q3_X_D - epsilon;
    T_lbin_val = t_eixo_CD_lbin;
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A) + R_By_lb * (x_in - Q3_X_B) + F_C_eng.Fy_lb * (x_in - Q3_X_C);
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A) + R_Bz_lb * (x_in - Q3_X_B) + F_C_eng.Fz_lb * (x_in - Q3_X_C);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

    // Ponto D (Roda Dentada)
    x_in = Q3_X_D;
    T_lbin_val = t_eixo_CD_lbin;
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A) + R_By_lb * (x_in - Q3_X_B) + F_C_eng.Fy_lb * (x_in - Q3_X_C);
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A) + R_Bz_lb * (x_in - Q3_X_B) + F_C_eng.Fz_lb * (x_in - Q3_X_C);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Antes
    T_lbin_val = t_eixo_DE_lbin;
    Vy_lbf_val += F_D_roda.Fy_lb;
    Vz_lbf_val += F_D_roda.Fz_lb;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Apos

    // Ponto antes de E
    x_in = Q3_X_E - epsilon;
    T_lbin_val = t_eixo_DE_lbin;
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A) + R_By_lb * (x_in - Q3_X_B) + F_C_eng.Fy_lb * (x_in - Q3_X_C) + F_D_roda.Fy_lb * (x_in - Q3_X_D);
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A) + R_Bz_lb * (x_in - Q3_X_B) + F_C_eng.Fz_lb * (x_in - Q3_X_C) + F_D_roda.Fz_lb * (x_in - Q3_X_D);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

    // Ponto E (Mancal)
    x_in = Q3_X_E;
    T_lbin_val = t_eixo_DE_lbin;                                                                                                                                                                                                     // Torque já deve ser zero ou próximo
    Mz_lbin_val = F_A.Fy_lb * (x_in - Q3_X_A) + R_By_lb * (x_in - Q3_X_B) + F_C_eng.Fy_lb * (x_in - Q3_X_C) + F_D_roda.Fy_lb * (x_in - Q3_X_D);                                                                                      // Deve ser ~0
    My_lbin_val = F_A.Fz_lb * (x_in - Q3_X_A) + R_Bz_lb * (x_in - Q3_X_B) + F_C_eng.Fz_lb * (x_in - Q3_X_C) + F_D_roda.Fz_lb * (x_in - Q3_X_D);                                                                                      // Deve ser ~0
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Antes
    Vy_lbf_val += R_Ey_lb;
    Vz_lbf_val += R_Ez_lb; // Cortantes devem zerar
    if (fabs(Vy_lbf_val) < epsilon)
        Mz_lbin_val = 0.0;
    if (fabs(Vz_lbf_val) < epsilon)
        My_lbin_val = 0.0;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON); // Apos

    // --- Impressao dos calculos de diametro em respostas.txt (Imperial com conversao) ---
    fprintf(fp_respostas, "\nd) Calculo dos Diametros Minimos do Eixo (em polegadas, com equivalentes SI):\n");
    double c_size_default = 0.85;
    double sn_linha_psi_default = q3_calcular_sn_linha(Q3_S_N_SURF_PSI, c_size_default, Q3_C_RELIAB_99);
    fprintf(fp_respostas, "   (Usando Sn_linha_padrao = %.0f psi (%.2f MPa), Sy = %.0f psi (%.2f MPa), N = %.1f)\n",
            sn_linha_psi_default, sn_linha_psi_default * PSI_PARA_PASCAL * PASCAL_PARA_MPA,
            Q3_S_Y_PSI, Q3_S_Y_PSI * PSI_PARA_PASCAL * PASCAL_PARA_MPA,
            Q3_FATOR_PROJETO_N);

    // Re-calculo dos momentos fletores resultantes para os pontos de interesse (A, B, C, D, E) em lb.in
    double m_res_A_lbin = 0.0; // No ponto da polia, momento é zero antes do primeiro mancal
    double d_A_in = q3_calcular_diametro_eixo_flexo_torcao(m_res_A_lbin, t_eixo_AB_lbin, Q3_KT_POLIA_A, sn_linha_psi_default, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto A (Polia, x=%.1f in): M_res=%.2f lb.in, T=%.2f lb.in, D_A = %.3f in (%.2f mm)\n", Q3_X_A, m_res_A_lbin, t_eixo_AB_lbin, d_A_in, d_A_in * POLEGADA_PARA_METRO * 1000);

    double mz_B_lbin = F_A.Fy_lb * (Q3_X_B - Q3_X_A);
    double my_B_lbin = F_A.Fz_lb * (Q3_X_B - Q3_X_A); // F_A.Fz_lb é 0
    double m_res_B_lbin = sqrt(pow(mz_B_lbin, 2) + pow(my_B_lbin, 2));
    double c_size_B = (d_A_in > 2.0 && d_A_in <= 3.0) ? 0.82 : 0.85; // Aproximacao
    double sn_linha_B_psi = q3_calcular_sn_linha(Q3_S_N_SURF_PSI, c_size_B, Q3_C_RELIAB_99);
    double d_B_in = q3_calcular_diametro_eixo_flexo_torcao(m_res_B_lbin, t_eixo_BC_lbin, Q3_KT_MANCAL_B, sn_linha_B_psi, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto B (Mancal, x=%.1f in): M_res=%.2f lb.in, T=%.2f lb.in, D_B = %.3f in (%.2f mm) (sn'=%.0f psi)\n", Q3_X_B, m_res_B_lbin, t_eixo_BC_lbin, d_B_in, d_B_in * POLEGADA_PARA_METRO * 1000, sn_linha_B_psi);

    double mz_C_lbin = F_A.Fy_lb * (Q3_X_C - Q3_X_A) + R_By_lb * (Q3_X_C - Q3_X_B);
    double my_C_lbin = F_A.Fz_lb * (Q3_X_C - Q3_X_A) + R_Bz_lb * (Q3_X_C - Q3_X_B); // F_A.Fz_lb é 0
    double m_res_C_lbin = sqrt(pow(mz_C_lbin, 2) + pow(my_C_lbin, 2));
    double d_C_in = q3_calcular_diametro_eixo_flexo_torcao(m_res_C_lbin, t_eixo_CD_lbin, Q3_KT_ENGRENAGEM_C, sn_linha_psi_default, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto C (Engrenagem, x=%.1f in): M_res=%.2f lb.in, T=%.2f lb.in, D_C = %.3f in (%.2f mm)\n", Q3_X_C, m_res_C_lbin, t_eixo_CD_lbin, d_C_in, d_C_in * POLEGADA_PARA_METRO * 1000);

    double mz_D_lbin = F_A.Fy_lb * (Q3_X_D - Q3_X_A) + R_By_lb * (Q3_X_D - Q3_X_B) + F_C_eng.Fy_lb * (Q3_X_D - Q3_X_C);
    double my_D_lbin = F_A.Fz_lb * (Q3_X_D - Q3_X_A) + R_Bz_lb * (Q3_X_D - Q3_X_B) + F_C_eng.Fz_lb * (Q3_X_D - Q3_X_C);
    double m_res_D_lbin = sqrt(pow(mz_D_lbin, 2) + pow(my_D_lbin, 2));
    double d_D_in = q3_calcular_diametro_eixo_flexo_torcao(m_res_D_lbin, t_eixo_DE_lbin, Q3_KT_RODA_DENTADA_D, sn_linha_psi_default, Q3_S_Y_PSI, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto D (Roda Dentada, x=%.1f in): M_res=%.2f lb.in, T=%.2f lb.in, D_D = %.3f in (%.2f mm)\n", Q3_X_D, m_res_D_lbin, t_eixo_DE_lbin, d_D_in, d_D_in * POLEGADA_PARA_METRO * 1000);

    double v_res_E_lbf = sqrt(pow(R_Ey_lb, 2) + pow(R_Ez_lb, 2));
    double d_E_in = q3_calcular_diametro_eixo_cisalhamento_puro(fabs(v_res_E_lbf), Q3_KT_MANCAL_E, sn_linha_psi_default, Q3_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto E (Mancal, x=%.1f in): V_res=%.2f lb, D_E = %.3f in (%.2f mm) (por cisalhamento)\n", Q3_X_E, fabs(v_res_E_lbf), d_E_in, d_E_in * POLEGADA_PARA_METRO * 1000);

    fprintf(fp_respostas, "\nDados para os diagramas (em SI) foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 3 ---\n\n\n");
}