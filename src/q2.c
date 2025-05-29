#include <stdio.h>
#include <math.h>
#include "q2.h" // Assumindo que q2.h declara void q2(FILE*, FILE*);

static double q2_para_radianos_static(double graus) {
    return graus * M_PI / 180.0;
}

// Função para calcular diâmetro (retorna em polegadas)
static double q2_calcular_diametro_eixo_ASME(double M_res_lbin, double T_lbin, double Kt,
                                           double sn_linha_psi, double sy_psi, double N_fator) {
    if (sn_linha_psi <= 0 || sy_psi <= 0) return -1.0; // Evita divisão por zero ou erro
    double termo_flexao_quad = pow((Kt * M_res_lbin / sn_linha_psi), 2.0);
    double termo_torcao_quad = (3.0 / 4.0) * pow((T_lbin / sy_psi), 2.0);

    if (fabs(M_res_lbin) < 1e-6 && fabs(T_lbin) < 1e-6) return 0.0; // Sem carga, diâmetro pode ser mínimo prático
    if (fabs(M_res_lbin) < 1e-6) termo_flexao_quad = 0; // Se M=0, Kt não se aplica a termo de flexão
    if (fabs(T_lbin) < 1e-6) termo_torcao_quad = 0;

    double dentro_raiz = sqrt(termo_flexao_quad + termo_torcao_quad);
    if (dentro_raiz < 1e-9 && (fabs(M_res_lbin) > 1e-6 || fabs(T_lbin) > 1e-6)) {
         // Caso raro onde M e T são pequenos mas não zero, e sn_linha ou sy são muito grandes
         // resultando em um valor quase zero dentro da raiz. Evita problemas com pow.
         // Retorna um diâmetro pequeno mas não zero.
         return 0.1; // Diâmetro mínimo prático arbitrário
    }
    double diametro_cubo = (32.0 * N_fator / M_PI) * dentro_raiz;
    if (diametro_cubo < 0) return -1.0; // Não deveria acontecer com sqrt
    return pow(diametro_cubo, 1.0 / 3.0);
}


static Q2_ForcaImperial q2_calcular_escrever_parametros_engrenagem_b(FILE *fp_respostas, double torque_lb_in) {
    Q2_ForcaImperial forca_B = {0.0, 0.0};
    fprintf(fp_respostas, "Engrenagem B (Questao 2 - recebe potencia, pinhão externo acima):\n");

    double d_gear_b_in = Q2_NUM_DENTES_GEARB / Q2_PASSO_DIAMETRAL_GEARB;
    fprintf(fp_respostas, "  Diametro Primitivo (D_gearB): %.2f in (%.4f m)\n",
            d_gear_b_in, d_gear_b_in * POLEGADA_PARA_METRO);

    double w_tb_lb = torque_lb_in / (d_gear_b_in / 2.0);
    forca_B.Fz_lb = w_tb_lb; 
    fprintf(fp_respostas, "  Forca Tangencial no eixo devido ao pinhão (F_Bz): %.2f lb (%.2f N) (Horizontal)\n",
            forca_B.Fz_lb, forca_B.Fz_lb * LBF_PARA_NEWTON);

    double angulo_pressao_rad = q2_para_radianos_static(Q2_ANGULO_PRESSAO_GEARB_DEG);
    double w_rb_lb = w_tb_lb * tan(angulo_pressao_rad);
    forca_B.Fy_lb = -w_rb_lb; 
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
    forca_D.Fy_lb = -f_polia_d_total_lb * cos(angulo_forca_rad); 
    forca_D.Fz_lb = f_polia_d_total_lb * sin(angulo_forca_rad);  
    fprintf(fp_respostas, "  Componente Vertical da Forca (F_Dy): %.2f lb (%.2f N) (Para Baixo)\n",
            forca_D.Fy_lb, forca_D.Fy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "  Componente Horizontal da Forca (F_Dz): %.2f lb (%.2f N) (Para Direita)\n",
            forca_D.Fz_lb, forca_D.Fz_lb * LBF_PARA_NEWTON);
    return forca_D;
}

void q2(FILE *fp_respostas, FILE *fp_dados) {
    fprintf(fp_respostas, "--- Questao 2 ---\n\n");
    // ... (impressão dos parâmetros da questão como antes) ...
    fprintf(fp_respostas, "Refazendo a Questao 1 com novos parametros:\n");
    fprintf(fp_respostas, "  Potencia: %.1f HP, Rotacao: %.0f RPM\n", Q2_POTENCIA_HP, Q2_ROTACAO_RPM);
    fprintf(fp_respostas, "  Engrenagem B: %d dentes, Passo Diametral %d\n", (int)Q2_NUM_DENTES_GEARB, (int)Q2_PASSO_DIAMETRAL_GEARB);
    fprintf(fp_respostas, "  Polia D: Diametro %.1f in (%.4f m)\n\n", Q2_DIAMETRO_POLIA_D_IN, Q2_DIAMETRO_POLIA_D_IN * POLEGADA_PARA_METRO);

    // --- a) Torque ---
    double torque_eixo_lb_in = (Q2_POTENCIA_HP * Q2_FATOR_TORQUE) / Q2_ROTACAO_RPM;
    // ... (impressão do torque em lb.in e N.m como antes) ...
    fprintf(fp_respostas, "a) Magnitude do Torque no Eixo (T_BD):\n");
    fprintf(fp_respostas, "  T = (%.1f hp * %.1f / %.0f rpm) = %.2f lb.in  (%.2f N.m)\n",
            Q2_POTENCIA_HP, Q2_FATOR_TORQUE, Q2_ROTACAO_RPM, torque_eixo_lb_in, torque_eixo_lb_in * LBIN_PARA_NM);
    fprintf(fp_respostas, "  Este torque atua entre a Engrenagem B (entrada) e a Polia D (saida).\n");
    fprintf(fp_respostas, "  Torque T_AB (antes de B) = 0 lb.in (0 N.m)\n");
    fprintf(fp_respostas, "  Torque T_DC (depois de D) = 0 lb.in (0 N.m)\n\n");

    // --- b) Forcas nos Elementos ---
    fprintf(fp_respostas, "b) Forcas nos Elementos de Transmissao:\n");
    Q2_ForcaImperial F_B = q2_calcular_escrever_parametros_engrenagem_b(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");
    Q2_ForcaImperial F_D = q2_calcular_escrever_parametros_polia_d(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");

    // --- c) Reacoes nos Mancais A e C ---
    fprintf(fp_respostas, "c) Reacoes nos Mancais (Mancal em A em x=%.1f in, Mancal em C em x=%.1f in):\n", Q2_X_A, Q2_X_C);
    double R_Ay_lb, R_Cy_lb, R_Az_lb, R_Cz_lb;
    // ... (cálculo das reações R_Ay_lb, etc. em imperial como antes) ...
    R_Cy_lb = -(F_B.Fy_lb * Q2_X_B + F_D.Fy_lb * Q2_X_D) / Q2_X_C;
    R_Ay_lb = -(F_B.Fy_lb + F_D.Fy_lb + R_Cy_lb);
    R_Cz_lb = -(F_B.Fz_lb * Q2_X_B + F_D.Fz_lb * Q2_X_D) / Q2_X_C;
    R_Az_lb = -(F_B.Fz_lb + F_D.Fz_lb + R_Cz_lb);
    // ... (impressão das reações em lb e N como antes) ...
    fprintf(fp_respostas, "  Plano Vertical (XY):\n");
    fprintf(fp_respostas, "    R_Ay = %.2f lb (%.2f N)\n", R_Ay_lb, R_Ay_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    R_Cy = %.2f lb (%.2f N)\n", R_Cy_lb, R_Cy_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "  Plano Horizontal (XZ):\n");
    fprintf(fp_respostas, "    R_Az = %.2f lb (%.2f N)\n", R_Az_lb, R_Az_lb * LBF_PARA_NEWTON);
    fprintf(fp_respostas, "    R_Cz = %.2f lb (%.2f N)\n\n", R_Cz_lb, R_Cz_lb * LBF_PARA_NEWTON);

    // --- d) Dados para Diagramas (Convertidos para SI ao escrever em dados.txt) ---
    // A lógica de geração de pontos e conversão para fp_dados permanece a mesma de q1.c,
    // usando as variáveis R_Ay_lb, F_B.Fy_lb, torque_eixo_lb_in, etc., da Questão 2.
    fprintf(fp_dados, "# Questao 2\n");
    // Formato: x[m] T[N.m] Mz[N.m] Vy[N] My[N.m] Vz[N]
    double x_in, T_lbin_val, Mz_lbin_val, Vy_lbf_val, My_lbin_val, Vz_lbf_val;
    double epsilon = 1e-6;

    // Ponto A (x=0)
    x_in = Q2_X_A; T_lbin_val = 0.0; Mz_lbin_val = 0.0; My_lbin_val = 0.0;
    Vy_lbf_val = 0.0; Vz_lbf_val = 0.0; 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    Vy_lbf_val = R_Ay_lb; Vz_lbf_val = R_Az_lb; 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto antes de B
    x_in = Q2_X_B - epsilon; T_lbin_val = 0.0;
    Mz_lbin_val = R_Ay_lb * x_in; My_lbin_val = R_Az_lb * x_in;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto B
    x_in = Q2_X_B;
    T_lbin_val = 0.0; Mz_lbin_val = R_Ay_lb * x_in; My_lbin_val = R_Az_lb * x_in; 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    T_lbin_val = torque_eixo_lb_in; 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    Vy_lbf_val += F_B.Fy_lb; Vz_lbf_val += F_B.Fz_lb; 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto antes de D
    x_in = Q2_X_D - epsilon; T_lbin_val = torque_eixo_lb_in;
    Mz_lbin_val = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q2_X_B);
    My_lbin_val = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q2_X_B);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto D
    x_in = Q2_X_D;
    T_lbin_val = torque_eixo_lb_in; Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q2_X_B); My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q2_X_B); 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    T_lbin_val = 0.0; 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);
    Vy_lbf_val += F_D.Fy_lb; Vz_lbf_val += F_D.Fz_lb; 
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto antes de C
    x_in = Q2_X_C - epsilon; T_lbin_val = 0.0;
    Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q2_X_B) + F_D.Fy_lb*(x_in-Q2_X_D);
    My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q2_X_B) + F_D.Fz_lb*(x_in-Q2_X_D);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON);

    // Ponto C (Mancal)
    x_in = Q2_X_C; T_lbin_val = 0.0;
    Mz_lbin_val = R_Ay_lb*x_in + F_B.Fy_lb*(x_in-Q2_X_B) + F_D.Fy_lb*(x_in-Q2_X_D);
    My_lbin_val = R_Az_lb*x_in + F_B.Fz_lb*(x_in-Q2_X_B) + F_D.Fz_lb*(x_in-Q2_X_D);
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON); 
    Vy_lbf_val += R_Cy_lb; Vz_lbf_val += R_Cz_lb; 
    if (fabs(Vy_lbf_val) < epsilon) Mz_lbin_val = 0.0;
    if (fabs(Vz_lbf_val) < epsilon) My_lbin_val = 0.0;
    fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in*POLEGADA_PARA_METRO, T_lbin_val*LBIN_PARA_NM, Mz_lbin_val*LBIN_PARA_NM, Vy_lbf_val*LBF_PARA_NEWTON, My_lbin_val*LBIN_PARA_NM, Vz_lbf_val*LBF_PARA_NEWTON); 

    // --- e) Cálculo dos Diâmetros Mínimos do Eixo (Aço SAE 1040 Estirado a Frio) ---
    fprintf(fp_respostas, "e) Calculo dos Diametros Minimos do Eixo (Aco SAE 1040 Estirado a Frio):\n");
    double sn_base_psi = Q2_S_N_BASE_PSI;
    double sn_linha_psi = sn_base_psi * Q2_C_SIZE_PADRAO * Q2_C_RELIAB_99;
    fprintf(fp_respostas, "   Material: Aco SAE 1040 Estirado a Frio\n");
    fprintf(fp_respostas, "     Su = %.0f psi (%.1f MPa)\n", Q2_S_U_PSI, Q2_S_U_PSI * PSI_PARA_PASCAL * PASCAL_PARA_MPA);
    fprintf(fp_respostas, "     Sy = %.0f psi (%.1f MPa)\n", Q2_S_Y_PSI, Q2_S_Y_PSI * PSI_PARA_PASCAL * PASCAL_PARA_MPA);
    fprintf(fp_respostas, "     Sn (base, usinado/est.frio) ~ %.0f psi (%.1f MPa)\n", sn_base_psi, sn_base_psi * PSI_PARA_PASCAL * PASCAL_PARA_MPA);
    fprintf(fp_respostas, "     Cs = %.2f (fator de tamanho padrao)\n", Q2_C_SIZE_PADRAO);
    fprintf(fp_respostas, "     CR = %.2f (para 99%% confiabilidade)\n", Q2_C_RELIAB_99);
    fprintf(fp_respostas, "     Sn' (corrigido) = %.0f psi (%.1f MPa)\n", sn_linha_psi, sn_linha_psi * PSI_PARA_PASCAL * PASCAL_PARA_MPA);
    fprintf(fp_respostas, "   Parametros de Projeto:\n");
    fprintf(fp_respostas, "     Fator de Seguranca (N) = %.1f\n", Q2_FATOR_PROJETO_N);

    // Pontos críticos para diâmetro: A, B, D, C
    // Momentos fletores resultantes (MR_lbin) e Torques (T_lbin_val_ponto) nesses pontos
    double MR_A_lbin, T_A_lbin_val, d_A_in;
    double MR_B_lbin, T_B_lbin_val, d_B_in;
    double MR_D_lbin, T_D_lbin_val, d_D_in;
    double MR_C_lbin, T_C_lbin_val, d_C_in;

    // Ponto A (Mancal, x=0)
    // Mz e My são 0 em A. Vy e Vz são as reações. M_res = 0.
    // Torque T_AB = 0
    MR_A_lbin = 0.0; T_A_lbin_val = 0.0; // Torque no segmento AB é 0
    d_A_in = q2_calcular_diametro_eixo_ASME(MR_A_lbin, T_A_lbin_val, Q2_KT_MANCAL, sn_linha_psi, Q2_S_Y_PSI, Q2_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto A (Mancal, x=%.1f in):\n", Q2_X_A);
    fprintf(fp_respostas, "    M_res = %.2f lb.in (%.2f N.m), T = %.2f lb.in (%.2f N.m), Kt = %.1f\n", MR_A_lbin, MR_A_lbin*LBIN_PARA_NM, T_A_lbin_val, T_A_lbin_val*LBIN_PARA_NM, Q2_KT_MANCAL);
    fprintf(fp_respostas, "    Diametro D_A = %.3f in (%.2f mm)\n", d_A_in, d_A_in * POLEGADA_PARA_METRO * 1000.0);

    // Ponto B (Engrenagem, x=Q2_X_B)
    Mz_lbin_val = R_Ay_lb * Q2_X_B; My_lbin_val = R_Az_lb * Q2_X_B;
    MR_B_lbin = sqrt(pow(Mz_lbin_val, 2) + pow(My_lbin_val, 2));
    T_B_lbin_val = torque_eixo_lb_in; // Torque T_BD
    d_B_in = q2_calcular_diametro_eixo_ASME(MR_B_lbin, T_B_lbin_val, Q2_KT_ENG_POLIA, sn_linha_psi, Q2_S_Y_PSI, Q2_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto B (Engrenagem, x=%.1f in):\n", Q2_X_B);
    fprintf(fp_respostas, "    M_res = %.2f lb.in (%.2f N.m), T = %.2f lb.in (%.2f N.m), Kt = %.1f\n", MR_B_lbin, MR_B_lbin*LBIN_PARA_NM, T_B_lbin_val, T_B_lbin_val*LBIN_PARA_NM, Q2_KT_ENG_POLIA);
    fprintf(fp_respostas, "    Diametro D_B = %.3f in (%.2f mm)\n", d_B_in, d_B_in * POLEGADA_PARA_METRO * 1000.0);

    // Ponto D (Polia, x=Q2_X_D)
    Mz_lbin_val = R_Ay_lb * Q2_X_D + F_B.Fy_lb * (Q2_X_D - Q2_X_B);
    My_lbin_val = R_Az_lb * Q2_X_D + F_B.Fz_lb * (Q2_X_D - Q2_X_B);
    MR_D_lbin = sqrt(pow(Mz_lbin_val, 2) + pow(My_lbin_val, 2));
    T_D_lbin_val = torque_eixo_lb_in; // Torque T_BD ainda atua em D antes da polia removê-lo
    d_D_in = q2_calcular_diametro_eixo_ASME(MR_D_lbin, T_D_lbin_val, Q2_KT_ENG_POLIA, sn_linha_psi, Q2_S_Y_PSI, Q2_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto D (Polia, x=%.1f in):\n", Q2_X_D);
    fprintf(fp_respostas, "    M_res = %.2f lb.in (%.2f N.m), T = %.2f lb.in (%.2f N.m), Kt = %.1f\n", MR_D_lbin, MR_D_lbin*LBIN_PARA_NM, T_D_lbin_val, T_D_lbin_val*LBIN_PARA_NM, Q2_KT_ENG_POLIA);
    fprintf(fp_respostas, "    Diametro D_D = %.3f in (%.2f mm)\n", d_D_in, d_D_in * POLEGADA_PARA_METRO * 1000.0);

    // Ponto C (Mancal, x=Q2_X_C)
    Mz_lbin_val = R_Ay_lb*Q2_X_C + F_B.Fy_lb*(Q2_X_C-Q2_X_B) + F_D.Fy_lb*(Q2_X_C-Q2_X_D); // Deveria ser ~0
    My_lbin_val = R_Az_lb*Q2_X_C + F_B.Fz_lb*(Q2_X_C-Q2_X_B) + F_D.Fz_lb*(Q2_X_C-Q2_X_D); // Deveria ser ~0
    if (fabs(Mz_lbin_val) < epsilon * 100) Mz_lbin_val = 0.0; // Ajuste para zerar se muito pequeno
    if (fabs(My_lbin_val) < epsilon * 100) My_lbin_val = 0.0;
    MR_C_lbin = sqrt(pow(Mz_lbin_val, 2) + pow(My_lbin_val, 2));
    T_C_lbin_val = 0.0; // Torque no segmento DC é 0
    d_C_in = q2_calcular_diametro_eixo_ASME(MR_C_lbin, T_C_lbin_val, Q2_KT_MANCAL, sn_linha_psi, Q2_S_Y_PSI, Q2_FATOR_PROJETO_N);
    fprintf(fp_respostas, "  Ponto C (Mancal, x=%.1f in):\n", Q2_X_C);
    fprintf(fp_respostas, "    M_res = %.2f lb.in (%.2f N.m), T = %.2f lb.in (%.2f N.m), Kt = %.1f\n", MR_C_lbin, MR_C_lbin*LBIN_PARA_NM, T_C_lbin_val, T_C_lbin_val*LBIN_PARA_NM, Q2_KT_MANCAL);
    fprintf(fp_respostas, "    Diametro D_C = %.3f in (%.2f mm)\n", d_C_in, d_C_in * POLEGADA_PARA_METRO * 1000.0);

    fprintf(fp_respostas, "\nDados para os diagramas (em SI) foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 2 ---\n\n\n");
}