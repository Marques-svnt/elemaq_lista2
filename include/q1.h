// Definicao de PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

void q1(FILE *fp_respostas, FILE *fp_dados);

// --- Fatores de Conversão ---
#define POLEGADA_PARA_METRO (0.0254)
#define LBF_PARA_NEWTON (4.4482216153)
#define LBIN_PARA_NM (LBF_PARA_NEWTON * POLEGADA_PARA_METRO) // ~0.1129848 N.m

// --- Constantes de entrada para a Questao 1 (Imperial) ---
#define Q1_POTENCIA_HP 30.0
#define Q1_ROTACAO_RPM 550.0
#define Q1_FATOR_TORQUE 63000.0

// Engrenagem B
#define Q1_PASSO_DIAMETRAL_GEARB 6.0
#define Q1_NUM_DENTES_GEARB 96.0
#define Q1_ANGULO_PRESSAO_GEARB_DEG 20.0

// Polia D
#define Q1_DIAMETRO_POLIA_D_IN 10.0
#define Q1_ANGULO_FORCA_POLIA_D_DEG 40.0 // Angulo com a vertical
#define Q1_FATOR_FORCA_POLIA_V 1.5

// Posicoes (polegadas) - Mancais em A e C
#define X_A 0.0
#define X_B 10.0 // Engrenagem B
#define X_D 20.0 // Polia D
#define X_C 26.0 // Mancal C

// Estrutura para armazenar forcas em um ponto (Fy, Fz em lb)
typedef struct
{
    double Fy_lb;
    double Fz_lb;
} Forca;
