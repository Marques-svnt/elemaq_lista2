#include <stdio.h>
#include <math.h>
#include "q2.h" // Incluir q2.h se houver declarações lá

// Funcao para converter graus em radianos
static double q2_para_radianos_static(double graus) { // Renomeada para evitar conflito se q1.c for compilado junto sem static
    return graus * M_PI / 180.0;
}

// Funcao para calcular e escrever os parametros da Engrenagem B para Q2
static Q2_Forca q2_calcular_escrever_parametros_engrenagem_b(FILE *fp_respostas, double torque_lb_in) {
    Q2_Forca forca_B = {0.0, 0.0};
    fprintf(fp_respostas, "Engrenagem B (Questao 2 - recebe potencia, pinhão externo acima):\n");

    double d_gear_b_in = Q2_NUM_DENTES_GEARB / Q2_PASSO_DIAMETRAL_GEARB;
    fprintf(fp_respostas, "  Diametro Primitivo (D_gearB): %.2f in\n", d_gear_b_in);

    double w_tb_lb = torque_lb_in / (d_gear_b_in / 2.0);
    forca_B.Fz = w_tb_lb; // Tangencial no eixo B (Horizontal)
    fprintf(fp_respostas, "  Forca Tangencial no eixo devido ao pinhão (F_Bz): %.2f lb (Horizontal)\n", forca_B.Fz);

    double angulo_pressao_rad = q2_para_radianos_static(Q2_ANGULO_PRESSAO_GEARB_DEG);
    double w_rb_lb = w_tb_lb * tan(angulo_pressao_rad);
    forca_B.Fy = -w_rb_lb; // Radial no eixo B (Vertical para Baixo)
    fprintf(fp_respostas, "  Forca Radial no eixo devido ao pinhão (F_By): %.2f lb (Vertical para Baixo)\n", forca_B.Fy);
    return forca_B;
}

// Funcao para calcular e escrever os parametros da Polia D para Q2
static Q2_Forca q2_calcular_escrever_parametros_polia_d(FILE *fp_respostas, double torque_lb_in) {
    Q2_Forca forca_D = {0.0, 0.0};
    fprintf(fp_respostas, "Polia D (Questao 2 - transmite potencia):\n");

    double f_efetiva_polia_d_lb = torque_lb_in / (Q2_DIAMETRO_POLIA_D_IN / 2.0);
    double f_polia_d_total_lb = Q2_FATOR_FORCA_POLIA_V * f_efetiva_polia_d_lb;
    fprintf(fp_respostas, "  Forca Resultante Total na Polia (F_poliaD): %.2f lb (estimada com fator %.1f)\n", f_polia_d_total_lb, Q2_FATOR_FORCA_POLIA_V);

    double angulo_forca_rad = q2_para_radianos_static(Q2_ANGULO_FORCA_POLIA_D_DEG);
    forca_D.Fy = -f_polia_d_total_lb * cos(angulo_forca_rad); // Componente vertical para baixo
    forca_D.Fz = f_polia_d_total_lb * sin(angulo_forca_rad);  // Componente horizontal para direita
    fprintf(fp_respostas, "  Componente Vertical da Forca (F_Dy): %.2f lb (Para Baixo)\n", forca_D.Fy);
    fprintf(fp_respostas, "  Componente Horizontal da Forca (F_Dz): %.2f lb (Para Direita)\n", forca_D.Fz);
    return forca_D;
}

void q2(FILE *fp_respostas, FILE *fp_dados) {
    fprintf(fp_respostas, "--- Questao 2 ---\n\n");
    fprintf(fp_respostas, "Refazendo a Questao 1 com novos parametros:\n");
    fprintf(fp_respostas, "  Potencia: %.1f HP, Rotacao: %.0f RPM\n", Q2_POTENCIA_HP, Q2_ROTACAO_RPM);
    fprintf(fp_respostas, "  Engrenagem B: %d dentes, Passo Diametral %d\n", (int)Q2_NUM_DENTES_GEARB, (int)Q2_PASSO_DIAMETRAL_GEARB);
    fprintf(fp_respostas, "  Polia D: Diametro %.1f in\n\n", Q2_DIAMETRO_POLIA_D_IN);

    // --- a) Torque ---
    double torque_eixo_lb_in = (Q2_POTENCIA_HP * Q2_FATOR_TORQUE) / Q2_ROTACAO_RPM;
    fprintf(fp_respostas, "a) Magnitude do Torque no Eixo (T_BD):\n");
    fprintf(fp_respostas, "  T = (%.1f * %.1f) / %.1f = %.2f lb.in\n",
            Q2_POTENCIA_HP, Q2_FATOR_TORQUE, Q2_ROTACAO_RPM, torque_eixo_lb_in);
    fprintf(fp_respostas, "  Este torque atua entre a Engrenagem B (entrada) e a Polia D (saida).\n");
    fprintf(fp_respostas, "  Torque T_AB (antes de B) = 0 lb.in\n");
    fprintf(fp_respostas, "  Torque T_DC (depois de D) = 0 lb.in\n\n");

    // --- b) Forcas nos Elementos ---
    fprintf(fp_respostas, "b) Forcas nos Elementos de Transmissao:\n");
    Q2_Forca F_B = q2_calcular_escrever_parametros_engrenagem_b(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");
    Q2_Forca F_D = q2_calcular_escrever_parametros_polia_d(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");

    // --- c) Reacoes nos Mancais A e C ---
    fprintf(fp_respostas, "c) Reacoes nos Mancais (Mancal em A em x=%.1f, Mancal em C em x=%.1f):\n", Q2_X_A, Q2_X_C);
    double R_Ay, R_Cy, R_Az, R_Cz;

    // Plano Vertical (XY)
    R_Cy = -(F_B.Fy * Q2_X_B + F_D.Fy * Q2_X_D) / Q2_X_C;
    R_Ay = -(F_B.Fy + F_D.Fy + R_Cy);
    fprintf(fp_respostas, "  Plano Vertical (XY):\n");
    fprintf(fp_respostas, "    R_Ay = %.2f lb\n", R_Ay);
    fprintf(fp_respostas, "    R_Cy = %.2f lb\n", R_Cy);

    // Plano Horizontal (XZ)
    R_Cz = -(F_B.Fz * Q2_X_B + F_D.Fz * Q2_X_D) / Q2_X_C;
    R_Az = -(F_B.Fz + F_D.Fz + R_Cz);
    fprintf(fp_respostas, "  Plano Horizontal (XZ):\n");
    fprintf(fp_respostas, "    R_Az = %.2f lb\n", R_Az);
    fprintf(fp_respostas, "    R_Cz = %.2f lb\n\n", R_Cz);

    // --- d) Dados para Diagramas (escritos em dados.txt) ---
    fprintf(fp_dados, "# Questao 2\n");
    // Formato: x T Mz Vy My Vz

    double x, T, Mz, Vy, My, Vz;
    double epsilon = 1e-6;

    // Ponto A (x=0)
    x = Q2_X_A;
    T = 0.0;
    Vy = R_Ay;
    Vz = R_Az;
    Mz = 0.0;
    My = 0.0;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, 0.0, My, 0.0); // Antes da reacao
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);  // Depois da reacao

    // Entre A e B (0 < x < X_B)
    x = Q2_X_B - epsilon;
    T = 0.0;
    Mz = R_Ay * x;
    My = R_Az * x;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);

    // Ponto B (x=X_B)
    x = Q2_X_B;
    Mz = R_Ay * x; // Mz e My antes do efeito de F_B para o momento no ponto B
    My = R_Az * x;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, 0.0, Mz, Vy, My, Vz); // Torque era 0 antes de B
    T = torque_eixo_lb_in; // Torque agora e aplicado
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);   // Torque agora e T_BD
    
    Vy += F_B.Fy; // Cortante apos F_B.Fy
    Vz += F_B.Fz; // Cortante apos F_B.Fz
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Cortante atualizada

    // Entre B e D (X_B < x < X_D)
    x = Q2_X_D - epsilon;
    T = torque_eixo_lb_in;
    Mz = R_Ay * x + F_B.Fy * (x - Q2_X_B);
    My = R_Az * x + F_B.Fz * (x - Q2_X_B);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);

    // Ponto D (x=X_D)
    x = Q2_X_D;
    Mz = R_Ay * x + F_B.Fy * (x - Q2_X_B); // Mz e My antes do efeito de F_D
    My = R_Az * x + F_B.Fz * (x - Q2_X_B);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Torque ainda e T_BD
    
    T = 0.0; // Torque agora e zero depois de D
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Torque agora e 0
    
    Vy += F_D.Fy; // Cortante apos F_D.Fy
    Vz += F_D.Fz; // Cortante apos F_D.Fz
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Cortante atualizada

    // Entre D e C (X_D < x < X_C)
    x = Q2_X_C - epsilon;
    T = 0.0;
    Mz = R_Ay * x + F_B.Fy * (x - Q2_X_B) + F_D.Fy * (x - Q2_X_D);
    My = R_Az * x + F_B.Fz * (x - Q2_X_B) + F_D.Fz * (x - Q2_X_D);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);

    // Ponto C (x=X_C) - Mancal
    x = Q2_X_C;
    T = 0.0;
    Mz = R_Ay * x + F_B.Fy * (x - Q2_X_B) + F_D.Fy * (x - Q2_X_D);
    My = R_Az * x + F_B.Fz * (x - Q2_X_B) + F_D.Fz * (x - Q2_X_D);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Antes da reacao R_C
    
    Vy += R_Cy; // Vy deve ir para perto de zero
    Vz += R_Cz; // Vz deve ir para perto de zero
    if (fabs(Vy) < epsilon) Mz = 0.0; // Forca Mz=0 se Vy~0
    if (fabs(Vz) < epsilon) My = 0.0; // Forca My=0 se Vz~0
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Depois da reacao R_C

    fprintf(fp_respostas, "Dados para os diagramas da Questao 2 foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 2 ---\n\n\n");
}