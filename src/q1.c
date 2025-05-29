#include <stdio.h>
#include <math.h>
#include "q1.h" // Assumindo que q1.h declara void q1(FILE*, FILE*); e typedef struct Forca;

static double q1_para_radianos_static(double graus)
{
        return graus * M_PI / 180.0;
}

static Q1_ForcaImperial q1_calcular_escrever_parametros_engrenagem_b(FILE *fp_respostas, double torque_lb_in)
{
        Q1_ForcaImperial forca_B = {0.0, 0.0};
        fprintf(fp_respostas, "Engrenagem B (recebe potencia, pinhão externo acima):\n");

        double d_gear_b_in = Q1_NUM_DENTES_GEARB / Q1_PASSO_DIAMETRAL_GEARB;
        fprintf(fp_respostas, "  Diametro Primitivo (D_gearB): %.2f in (%.4f m)\n",
                d_gear_b_in, d_gear_b_in * POLEGADA_PARA_METRO);

        double w_tb_lb = torque_lb_in / (d_gear_b_in / 2.0);
        forca_B.Fz_lb = w_tb_lb; // Horizontal
        fprintf(fp_respostas, "  Forca Tangencial no eixo (F_Bz): %.2f lb (%.2f N) (Horizontal)\n",
                forca_B.Fz_lb, forca_B.Fz_lb * LBF_PARA_NEWTON);

        double angulo_pressao_rad = q1_para_radianos_static(Q1_ANGULO_PRESSAO_GEARB_DEG);
        double w_rb_lb = w_tb_lb * tan(angulo_pressao_rad);
        forca_B.Fy_lb = -w_rb_lb; // Vertical para Baixo
        fprintf(fp_respostas, "  Forca Radial no eixo (F_By): %.2f lb (%.2f N) (Vertical para Baixo)\n",
                forca_B.Fy_lb, forca_B.Fy_lb * LBF_PARA_NEWTON);
        return forca_B;
}

static Q1_ForcaImperial q1_calcular_escrever_parametros_polia_d(FILE *fp_respostas, double torque_lb_in)
{
        Q1_ForcaImperial forca_D = {0.0, 0.0};
        fprintf(fp_respostas, "Polia D (transmite potencia):\n");
        fprintf(fp_respostas, "  Diametro da Polia D: %.2f in (%.4f m)\n",
                Q1_DIAMETRO_POLIA_D_IN, Q1_DIAMETRO_POLIA_D_IN * POLEGADA_PARA_METRO);

        double f_efetiva_polia_d_lb = torque_lb_in / (Q1_DIAMETRO_POLIA_D_IN / 2.0);
        double f_polia_d_total_lb = Q1_FATOR_FORCA_POLIA_V * f_efetiva_polia_d_lb;
        fprintf(fp_respostas, "  Forca Resultante Total na Polia (F_poliaD): %.2f lb (%.2f N) (estimada com fator %.1f)\n",
                f_polia_d_total_lb, f_polia_d_total_lb * LBF_PARA_NEWTON, Q1_FATOR_FORCA_POLIA_V);

        double angulo_forca_rad = q1_para_radianos_static(Q1_ANGULO_FORCA_POLIA_D_DEG);
        forca_D.Fy_lb = -f_polia_d_total_lb * cos(angulo_forca_rad); // Para Baixo
        forca_D.Fz_lb = f_polia_d_total_lb * sin(angulo_forca_rad);  // Para Direita
        fprintf(fp_respostas, "  Componente Vertical da Forca (F_Dy): %.2f lb (%.2f N) (Para Baixo)\n",
                forca_D.Fy_lb, forca_D.Fy_lb * LBF_PARA_NEWTON);
        fprintf(fp_respostas, "  Componente Horizontal da Forca (F_Dz): %.2f lb (%.2f N) (Para Direita)\n",
                forca_D.Fz_lb, forca_D.Fz_lb * LBF_PARA_NEWTON);
        return forca_D;
}

// Funcoes de calculo de diametro (Imperial)
static double q1_calcular_sn_linha(double sn_base_psi, double c_size, double c_reliab)
{
        return sn_base_psi * c_size * c_reliab;
}
static double q1_calcular_diametro_eixo_flexo_torcao(double M_res_lbin, double T_lbin, double Kt,
                                                     double sn_linha_psi, double sy_psi, double N_fator)
{
        if (sn_linha_psi <= 0 || sy_psi <= 0)
                return -1.0;
        double tf_sq = pow((Kt * M_res_lbin / sn_linha_psi), 2.0);
        double tt_sq = (3.0 / 4.0) * pow((T_lbin / sy_psi), 2.0);
        if (fabs(M_res_lbin) < 1e-6 && fabs(T_lbin) < 1e-6)
                return 0.0; // Sem carga
        if (fabs(M_res_lbin) < 1e-6)
                tf_sq = 0; // Apenas torção (simplificado)
        if (fabs(T_lbin) < 1e-6)
                tt_sq = 0; // Apenas flexão
        double dentro_raiz = sqrt(tf_sq + tt_sq);
        double diametro_cubo = (32.0 * N_fator / M_PI) * dentro_raiz;
        return pow(diametro_cubo, 1.0 / 3.0); // Retorna diâmetro em polegadas
}

void q1(FILE *fp_respostas, FILE *fp_dados)
{
        fprintf(fp_respostas, "--- Questao 1 ---\n\n");
        fprintf(fp_respostas, "Material do Eixo: Aco SAE 1040 Estirado a Frio\n");
        fprintf(fp_respostas, "a) Magnitude do torque no eixo em todos os pontos;\n");
        fprintf(fp_respostas, "b) Forcas nos elementos transmissores de potencia;\n");
        fprintf(fp_respostas, "c) Reacoes nos rolamentos;\n");
        fprintf(fp_respostas, "d) Diagramas de carga, cisalhamento, momento fletor e calculo de diametros.\n\n");

        // --- a) Torque (Calculado em Imperial) ---
        double torque_eixo_lb_in = (Q1_POTENCIA_HP * Q1_FATOR_TORQUE) / Q1_ROTACAO_RPM;
        fprintf(fp_respostas, "Torque no Eixo (T_BD):\n");
        fprintf(fp_respostas, "  T = (%.1f hp * %.1f / %.1f rpm) = %.2f lb.in  (%.2f N.m)\n",
                Q1_POTENCIA_HP, Q1_FATOR_TORQUE, Q1_ROTACAO_RPM, torque_eixo_lb_in, torque_eixo_lb_in * LBIN_PARA_NM);
        fprintf(fp_respostas, "  Este torque atua entre a Engrenagem B (entrada) e a Polia D (saida).\n");
        fprintf(fp_respostas, "  Torque T_AB (antes de B) = 0 lb.in (0 N.m)\n");
        fprintf(fp_respostas, "  Torque T_DC (depois de D) = 0 lb.in (0 N.m)\n\n");

        // --- b) Forcas nos Elementos (Calculadas em Imperial, impressas com conversao) ---
        fprintf(fp_respostas, "Forcas nos Elementos de Transmissao:\n");
        Q1_ForcaImperial F_B = q1_calcular_escrever_parametros_engrenagem_b(fp_respostas, torque_eixo_lb_in);
        fprintf(fp_respostas, "\n");
        Q1_ForcaImperial F_D = q1_calcular_escrever_parametros_polia_d(fp_respostas, torque_eixo_lb_in);
        fprintf(fp_respostas, "\n");

        // --- c) Reacoes nos Mancais A e C (Calculadas em Imperial, impressas com conversao) ---
        fprintf(fp_respostas, "Reacoes nos Mancais (Mancal em A em x=%.1f in, Mancal em C em x=%.1f in):\n", Q1_X_A, Q1_X_C);
        double R_Ay_lb, R_Cy_lb, R_Az_lb, R_Cz_lb;

        R_Cy_lb = -(F_B.Fy_lb * Q1_X_B + F_D.Fy_lb * Q1_X_D) / Q1_X_C;
        R_Ay_lb = -(F_B.Fy_lb + F_D.Fy_lb + R_Cy_lb);
        fprintf(fp_respostas, "  Plano Vertical (XY):\n");
        fprintf(fp_respostas, "    R_Ay = %.2f lb (%.2f N)\n", R_Ay_lb, R_Ay_lb * LBF_PARA_NEWTON);
        fprintf(fp_respostas, "    R_Cy = %.2f lb (%.2f N)\n", R_Cy_lb, R_Cy_lb * LBF_PARA_NEWTON);

        R_Cz_lb = -(F_B.Fz_lb * Q1_X_B + F_D.Fz_lb * Q1_X_D) / Q1_X_C;
        R_Az_lb = -(F_B.Fz_lb + F_D.Fz_lb + R_Cz_lb);
        fprintf(fp_respostas, "  Plano Horizontal (XZ):\n");
        fprintf(fp_respostas, "    R_Az = %.2f lb (%.2f N)\n", R_Az_lb, R_Az_lb * LBF_PARA_NEWTON);
        fprintf(fp_respostas, "    R_Cz = %.2f lb (%.2f N)\n\n", R_Cz_lb, R_Cz_lb * LBF_PARA_NEWTON);

        // --- d) Calculo de Diametros (Imperial, impressos com conversao) ---
        fprintf(fp_respostas, "Calculo dos Diametros Minimos do Eixo (Aco SAE 1040 Estirado a Frio):\n");
        double sn_linha_psi_default = q1_calcular_sn_linha(Q1_S_N_BASE_PSI, Q1_C_SIZE_DEFAULT, Q1_C_RELIAB_99);
        fprintf(fp_respostas, "   (Usando Sn_linha = %.0f psi (%.2f MPa), Sy = %.0f psi (%.2f MPa), N = %.1f)\n",
                sn_linha_psi_default, sn_linha_psi_default * PSI_PARA_PASCAL * PASCAL_PARA_MPA,
                Q1_S_Y_PSI, Q1_S_Y_PSI * PSI_PARA_PASCAL * PASCAL_PARA_MPA,
                Q1_FATOR_PROJETO_N);

        // Momentos fletores resultantes nos pontos B e D (em lb.in)
        double Mz_B_lbin = R_Ay_lb * Q1_X_B;
        double My_B_lbin = R_Az_lb * Q1_X_B;
        double M_res_B_lbin = sqrt(pow(Mz_B_lbin, 2) + pow(My_B_lbin, 2));

        double Mz_D_lbin = R_Ay_lb * Q1_X_D + F_B.Fy_lb * (Q1_X_D - Q1_X_B);
        double My_D_lbin = R_Az_lb * Q1_X_D + F_B.Fz_lb * (Q1_X_D - Q1_X_B);
        double M_res_D_lbin = sqrt(pow(Mz_D_lbin, 2) + pow(My_D_lbin, 2));

        // Diâmetro no ponto B (sob engrenagem, com torque T_BD, Kt para chaveta)
        double d_B_in = q1_calcular_diametro_eixo_flexo_torcao(M_res_B_lbin, torque_eixo_lb_in, Q1_KT_CHAVETA, sn_linha_psi_default, Q1_S_Y_PSI, Q1_FATOR_PROJETO_N);
        fprintf(fp_respostas, "  Ponto B (Engrenagem, x=%.1f in): M_res=%.2f lb.in, T=%.2f lb.in, Kt=%.1f, D_B = %.3f in (%.2f mm)\n",
                Q1_X_B, M_res_B_lbin, torque_eixo_lb_in, Q1_KT_CHAVETA, d_B_in, d_B_in * POLEGADA_PARA_METRO * 1000);

        // Diâmetro no ponto D (sob polia, com torque T_BD, Kt para chaveta)
        // O torque em D é o mesmo que em B para este problema (T_BD)
        double d_D_in = q1_calcular_diametro_eixo_flexo_torcao(M_res_D_lbin, torque_eixo_lb_in, Q1_KT_CHAVETA, sn_linha_psi_default, Q1_S_Y_PSI, Q1_FATOR_PROJETO_N);
        fprintf(fp_respostas, "  Ponto D (Polia, x=%.1f in): M_res=%.2f lb.in, T=%.2f lb.in, Kt=%.1f, D_D = %.3f in (%.2f mm)\n",
                Q1_X_D, M_res_D_lbin, torque_eixo_lb_in, Q1_KT_CHAVETA, d_D_in, d_D_in * POLEGADA_PARA_METRO * 1000);

        // (Opcional) Diâmetros nos mancais A e C (considerando apenas flexão ou um Kt para filete)
        // Em A (x=0), M=0, T=0 (antes de B). Diâmetro seria mínimo ou prático.
        // Em C (x=26), M=0 (deveria ser), T=0. Diâmetro seria mínimo ou prático.
        // Para simplificar, não calcularemos diâmetros em A e C aqui, pois os pontos críticos são B e D.
        fprintf(fp_respostas, "\n");

        // --- Dados para Diagramas (Convertidos para SI ao escrever em dados.txt) ---
        fprintf(fp_dados, "# Questao 1\n");
        // Formato: x[m] T[N.m] Mz[N.m] Vy[N] My[N.m] Vz[N]
        // (A lógica de escrita para dados.txt permanece a mesma da sua versão anterior,
        //  pois ela já converte os valores imperiais para SI antes de fprintf)

        double x_in, T_lbin_val, Mz_lbin_val, Vy_lbf_val, My_lbin_val, Vz_lbf_val;
        double epsilon = 1e-6;

        // Ponto A (x=0)
        x_in = Q1_X_A;
        T_lbin_val = 0.0;
        Mz_lbin_val = 0.0;
        My_lbin_val = 0.0;
        Vy_lbf_val = 0.0;
        Vz_lbf_val = 0.0;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);
        Vy_lbf_val = R_Ay_lb;
        Vz_lbf_val = R_Az_lb;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

        // Ponto antes de B
        x_in = Q1_X_B - epsilon;
        T_lbin_val = 0.0;
        Mz_lbin_val = R_Ay_lb * x_in;
        My_lbin_val = R_Az_lb * x_in;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

        // Ponto B
        x_in = Q1_X_B;
        T_lbin_val = 0.0;
        Mz_lbin_val = R_Ay_lb * x_in;
        My_lbin_val = R_Az_lb * x_in;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);
        T_lbin_val = torque_eixo_lb_in;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);
        Vy_lbf_val += F_B.Fy_lb;
        Vz_lbf_val += F_B.Fz_lb;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

        // Ponto antes de D
        x_in = Q1_X_D - epsilon;
        T_lbin_val = torque_eixo_lb_in;
        Mz_lbin_val = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q1_X_B);
        My_lbin_val = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q1_X_B);
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

        // Ponto D
        x_in = Q1_X_D;
        T_lbin_val = torque_eixo_lb_in;
        Mz_lbin_val = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q1_X_B);
        My_lbin_val = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q1_X_B);
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);
        T_lbin_val = 0.0;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);
        Vy_lbf_val += F_D.Fy_lb;
        Vz_lbf_val += F_D.Fz_lb;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

        // Ponto antes de C (Mancal)
        x_in = Q1_X_C - epsilon;
        T_lbin_val = 0.0;
        Mz_lbin_val = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q1_X_B) + F_D.Fy_lb * (x_in - Q1_X_D);
        My_lbin_val = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q1_X_B) + F_D.Fz_lb * (x_in - Q1_X_D);
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

        // Ponto C (Mancal)
        x_in = Q1_X_C;
        T_lbin_val = 0.0;
        Mz_lbin_val = R_Ay_lb * x_in + F_B.Fy_lb * (x_in - Q1_X_B) + F_D.Fy_lb * (x_in - Q1_X_D);
        My_lbin_val = R_Az_lb * x_in + F_B.Fz_lb * (x_in - Q1_X_B) + F_D.Fz_lb * (x_in - Q1_X_D);
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);
        Vy_lbf_val += R_Cy_lb;
        Vz_lbf_val += R_Cz_lb;
        if (fabs(Vy_lbf_val) < epsilon)
                Mz_lbin_val = 0.0;
        if (fabs(Vz_lbf_val) < epsilon)
                My_lbin_val = 0.0;
        fprintf(fp_dados, "%.4f %.2f %.2f %.2f %.2f %.2f\n", x_in * POLEGADA_PARA_METRO, T_lbin_val * LBIN_PARA_NM, Mz_lbin_val * LBIN_PARA_NM, Vy_lbf_val * LBF_PARA_NEWTON, My_lbin_val * LBIN_PARA_NM, Vz_lbf_val * LBF_PARA_NEWTON);

        fprintf(fp_respostas, "Dados para os diagramas (em SI) foram escritos em dados.txt.\n");
        fprintf(fp_respostas, "--- Fim da Questao 1 ---\n\n\n");
}