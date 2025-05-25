#include <stdio.h>
#include <math.h>
#include "q1.h" // Incluir q1.h se houver declarações lá

// Funcao para converter graus em radianos
static double q1_para_radianos(double graus) {
    return graus * M_PI / 180.0;
}

// Funcao para calcular e escrever os parametros da Engrenagem B
static Forca calcular_escrever_parametros_engrenagem_b(FILE *fp_respostas, double torque_lb_in) {
    Forca forca_B = {0.0, 0.0};
    fprintf(fp_respostas, "Engrenagem B (recebe potencia, pinhão externo acima):\n");

    double d_gear_b_in = Q1_NUM_DENTES_GEARB / Q1_PASSO_DIAMETRAL_GEARB;
    fprintf(fp_respostas, "  Diametro Primitivo (D_gearB): %.2f in\n", d_gear_b_in);

    double w_tb_lb = torque_lb_in / (d_gear_b_in / 2.0);
    // Forca tangencial exercida PELO PINHAO EXTERNO EM B (no eixo) e horizontal
    forca_B.Fz = w_tb_lb; // Assumindo +Z (precisaria definir sentido da rotacao para ser exato)
    fprintf(fp_respostas, "  Forca Tangencial no eixo devido ao pinhão (F_Bz): %.2f lb (Horizontal)\n", forca_B.Fz);

    double angulo_pressao_rad = q1_para_radianos(Q1_ANGULO_PRESSAO_GEARB_DEG);
    double w_rb_lb = w_tb_lb * tan(angulo_pressao_rad);
    // Forca radial exercida PELO PINHAO EXTERNO EM B (no eixo) e vertical para baixo
    forca_B.Fy = -w_rb_lb;
    fprintf(fp_respostas, "  Forca Radial no eixo devido ao pinhão (F_By): %.2f lb (Vertical para Baixo)\n", forca_B.Fy);
    return forca_B;
}

// Funcao para calcular e escrever os parametros da Polia D
static Forca calcular_escrever_parametros_polia_d(FILE *fp_respostas, double torque_lb_in) {
    Forca forca_D = {0.0, 0.0};
    fprintf(fp_respostas, "Polia D (transmite potencia):\n");

    double f_efetiva_polia_d_lb = torque_lb_in / (Q1_DIAMETRO_POLIA_D_IN / 2.0);
    double f_polia_d_total_lb = Q1_FATOR_FORCA_POLIA_V * f_efetiva_polia_d_lb;
    fprintf(fp_respostas, "  Forca Resultante Total na Polia (F_poliaD): %.2f lb (estimada com fator %.1f)\n", f_polia_d_total_lb, Q1_FATOR_FORCA_POLIA_V);

    // Forca a 40 graus com a vertical, para baixo e para a direita
    double angulo_forca_rad = q1_para_radianos(Q1_ANGULO_FORCA_POLIA_D_DEG);
    forca_D.Fy = -f_polia_d_total_lb * cos(angulo_forca_rad); // Componente vertical para baixo
    forca_D.Fz = f_polia_d_total_lb * sin(angulo_forca_rad);  // Componente horizontal para direita
    fprintf(fp_respostas, "  Componente Vertical da Forca (F_Dy): %.2f lb (Para Baixo)\n", forca_D.Fy);
    fprintf(fp_respostas, "  Componente Horizontal da Forca (F_Dz): %.2f lb (Para Direita)\n", forca_D.Fz);
    return forca_D;
}

void q1(FILE *fp_respostas, FILE *fp_dados) {
    fprintf(fp_respostas, "--- Questao 1 ---\n\n");
    fprintf(fp_respostas, "a) Magnitude do torque no eixo em todos os pontos;\n");
    fprintf(fp_respostas, "b) Calcular as forcas que atuam no eixo em todos os elementos transmissores de potencia;\n");
    fprintf(fp_respostas, "c) Calcular as reacoes nos rolamentos;\n");
    fprintf(fp_respostas, "d) Desenhar os diagramas completos de carga, cisalhamento e momento fletor.\n\n");

    // --- a) Torque ---
    double torque_eixo_lb_in = (Q1_POTENCIA_HP * Q1_FATOR_TORQUE) / Q1_ROTACAO_RPM;
    fprintf(fp_respostas, "Torque no Eixo (T_BD):\n");
    fprintf(fp_respostas, "  T = (P_hp * 63000) / n_rpm = (%.1f * %.1f) / %.1f = %.2f lb.in\n",
            Q1_POTENCIA_HP, Q1_FATOR_TORQUE, Q1_ROTACAO_RPM, torque_eixo_lb_in);
    fprintf(fp_respostas, "  Este torque atua entre a Engrenagem B (entrada) e a Polia D (saida).\n");
    fprintf(fp_respostas, "  Torque T_AB (antes de B) = 0 lb.in\n");
    fprintf(fp_respostas, "  Torque T_DC (depois de D) = 0 lb.in\n\n");

    // --- b) Forcas nos Elementos ---
    fprintf(fp_respostas, "Forcas nos Elementos de Transmissao:\n");
    Forca F_B = calcular_escrever_parametros_engrenagem_b(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");
    Forca F_D = calcular_escrever_parametros_polia_d(fp_respostas, torque_eixo_lb_in);
    fprintf(fp_respostas, "\n");

    // --- c) Reacoes nos Mancais A e C ---
    fprintf(fp_respostas, "Reacoes nos Mancais (Mancal em A em x=%.1f, Mancal em C em x=%.1f):\n", X_A, X_C);
    double R_Ay, R_Cy, R_Az, R_Cz;

    // Plano Vertical (XY) - Forcas: F_B.Fy em X_B, F_D.Fy em X_D
    // Sum M_A_z = 0: (F_B.Fy * X_B) + (F_D.Fy * X_D) + (R_Cy * X_C) = 0
    R_Cy = -(F_B.Fy * X_B + F_D.Fy * X_D) / X_C;
    // Sum F_y = 0: R_Ay + F_B.Fy + F_D.Fy + R_Cy = 0
    R_Ay = -(F_B.Fy + F_D.Fy + R_Cy);

    fprintf(fp_respostas, "  Plano Vertical (XY):\n");
    fprintf(fp_respostas, "    R_Ay = %.2f lb\n", R_Ay);
    fprintf(fp_respostas, "    R_Cy = %.2f lb\n", R_Cy);

    // Plano Horizontal (XZ) - Forcas: F_B.Fz em X_B, F_D.Fz em X_D
    // Sum M_A_y = 0: (F_B.Fz * X_B) + (F_D.Fz * X_D) + (R_Cz * X_C) = 0
    R_Cz = -(F_B.Fz * X_B + F_D.Fz * X_D) / X_C;
    // Sum F_z = 0: R_Az + F_B.Fz + F_D.Fz + R_Cz = 0
    R_Az = -(F_B.Fz + F_D.Fz + R_Cz);

    fprintf(fp_respostas, "  Plano Horizontal (XZ):\n");
    fprintf(fp_respostas, "    R_Az = %.2f lb\n", R_Az);
    fprintf(fp_respostas, "    R_Cz = %.2f lb\n\n", R_Cz);

    // --- d) Dados para Diagramas (escritos em dados.txt) ---
    fprintf(fp_dados, "# Questao 1\n");
    // Formato: x T Mz Vy My Vz

    double x, T, Mz, Vy, My, Vz;
    double epsilon = 1e-6; // Pequeno valor para pontos "antes" e "depois"

    // Ponto A (x=0)
    x = X_A;
    T = 0.0;
    Vy = R_Ay;        // Cortante logo apos o mancal A
    Vz = R_Az;        // Cortante logo apos o mancal A
    Mz = 0.0;
    My = 0.0;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, 0.0, My, 0.0); // Antes da reacao
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);  // Depois da reacao

    // Entre A e B (0 < x < 10)
    x = X_B - epsilon; // Ponto infinitesimalmente antes de B
    T = 0.0;
    // Vy e Vz constantes = R_Ay, R_Az
    Mz = R_Ay * x;
    My = R_Az * x;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);

    // Ponto B (x=10) - Entrada de Torque e Forcas F_B
    x = X_B;
    // Antes da Forca F_B (mas depois da entrada de torque)
    T = torque_eixo_lb_in; // Torque agora e aplicado
    // Mz, My continuam do ponto anterior
    Mz = R_Ay * x;
    My = R_Az * x;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, 0.0, Mz, Vy, My, Vz); // Torque era 0 antes de B
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);   // Torque agora e T_BD
    
    // Depois da Forca F_B
    Vy += F_B.Fy;
    Vz += F_B.Fz;
    // Mz e My nao mudam instantaneamente no ponto da forca para o diagrama de momento
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);

    // Entre B e D (10 < x < 20)
    x = X_D - epsilon; // Ponto infinitesimalmente antes de D
    T = torque_eixo_lb_in;
    // Vy e Vz constantes = R_Ay + F_B.Fy, R_Az + F_B.Fz
    Mz = R_Ay * x + F_B.Fy * (x - X_B);
    My = R_Az * x + F_B.Fz * (x - X_B);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);

    // Ponto D (x=20) - Saida de Torque e Forcas F_D
    x = X_D;
    // Antes da Forca F_D (mas antes da saida de torque)
    T = torque_eixo_lb_in;
    Mz = R_Ay * x + F_B.Fy * (x - X_B);
    My = R_Az * x + F_B.Fz * (x - X_B);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Torque ainda e T_BD
    
    T = 0.0; // Torque agora e zero depois de D
    // Depois da Forca F_D
    Vy += F_D.Fy;
    Vz += F_D.Fz;
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Torque agora e 0

    // Entre D e C (20 < x < 26)
    x = X_C - epsilon; // Ponto infinitesimalmente antes de C
    T = 0.0;
    // Vy e Vz constantes
    Mz = R_Ay * x + F_B.Fy * (x - X_B) + F_D.Fy * (x - X_D);
    My = R_Az * x + F_B.Fz * (x - X_B) + F_D.Fz * (x - X_D);
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz);

    // Ponto C (x=26) - Mancal
    x = X_C;
    T = 0.0;
    Mz = R_Ay * x + F_B.Fy * (x - X_B) + F_D.Fy * (x - X_D); // Momento deve ser proximo de zero se R_Cy balancear
    My = R_Az * x + F_B.Fz * (x - X_B) + F_D.Fz * (x - X_D); // Momento deve ser proximo de zero se R_Cz balancear
    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Antes da reacao R_C
    
    Vy += R_Cy; // Vy deve ir para perto de zero
    Vz += R_Cz; // Vz deve ir para perto de zero
    // Para fins de diagrama, o momento fletor no mancal que zera a cortante também zera o momento.
    // A maneira como Mz e My foram calculados já deve resultar em valores próximos de 0 em x=X_C
    // se as reações estiverem corretas. Forçar Mz=0 e My=0 aqui se Vy e Vz zerarem.
    if (fabs(Vy) < epsilon) Mz = 0.0;
    if (fabs(Vz) < epsilon) My = 0.0;

    fprintf(fp_dados, "%.2f %.2f %.2f %.2f %.2f %.2f\n", x, T, Mz, Vy, My, Vz); // Depois da reacao R_C

    fprintf(fp_respostas, "Dados para os diagramas foram escritos em dados.txt.\n");
    fprintf(fp_respostas, "--- Fim da Questao 1 ---\n\n\n");
}