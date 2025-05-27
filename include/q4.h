#ifndef Q4_H
#define Q4_H

void q4(FILE *fp_respostas, FILE *fp_dados);
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// --- Fatores de Conversão ---
#define POLEGADA_PARA_METRO (0.0254)
#define LBF_PARA_NEWTON (4.4482216153)
#define LBIN_PARA_NM (LBF_PARA_NEWTON * POLEGADA_PARA_METRO)
#define PSI_PARA_PASCAL (LBF_PARA_NEWTON / (POLEGADA_PARA_METRO * POLEGADA_PARA_METRO))
#define PASCAL_PARA_MPA (1.0 / 1000000.0)

// --- Constantes da Questao 4 (Imperial) ---
#define Q4_VELOCIDADE_RPM 480.0
#define Q4_FATOR_TORQUE 63000.0 // lb.in * rpm / hp

// Potencias (hp)
#define Q4_POTENCIA_CORR_C_HP 11.0 // Entrada
#define Q4_POTENCIA_ENG_B_HP 5.0   // Saida
#define Q4_POTENCIA_POLIA_D_HP 3.0 // Saida
#define Q4_POTENCIA_POLIA_E_HP 3.0 // Saida

// Engrenagem B
#define Q4_DIAMETRO_ENG_B_IN 3.0
#define Q4_ANGULO_PRESSAO_ENG_B_DEG 20.0

// Corrente Dentada C
#define Q4_DIAMETRO_CORR_C_IN 10.0
#define Q4_ANGULO_FORCA_CORR_C_DEG 15.0 // Com a VERTICAL, para baixo e esquerda (conforme Figura 3)

// Polia D
#define Q4_DIAMETRO_POLIA_D_IN 4.0
#define Q4_ANGULO_FORCA_POLIA_D_DEG 30.0 // Com a VERTICAL, para baixo e direita (conforme Figura 3)

// Polia E
#define Q4_DIAMETRO_POLIA_E_IN 4.0
// Forca horizontal para a direita

// Fator para forca em correia V
#define Q4_FATOR_FORCA_POLIA_V 1.5

// Posicoes Axiais dos Elementos (polegadas)
#define Q4_X_A 0.0  // Mancal A
#define Q4_X_B 5.0  // Engrenagem B
#define Q4_X_C 13.0 // Corrente C
#define Q4_X_D 20.0 // Polia D
#define Q4_X_E 26.0 // Polia E
#define Q4_X_F 31.0 // Mancal F

// Propriedades do Material (SAE 1137 OQT 1300) e Projeto (Imperial)
#define Q4_S_U_PSI 105000.0
#define Q4_S_Y_PSI 80000.0
#define Q4_S_N_BASE_PSI 45000.0
#define Q4_C_SIZE 0.80
#define Q4_C_RELIAB_99 0.81
#define Q4_FATOR_PROJETO_N 2.5
#define Q4_KT_ANEL_RETENCAO 3.0

// Estrutura para forcas (em lb)
typedef struct {
    double Fy_lb;
    double Fz_lb;
} Q4_ForcaImperial;



#endif 