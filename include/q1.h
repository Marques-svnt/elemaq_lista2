// Definicao de PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

void q1(FILE *fp_respostas, FILE *fp_dados);
// Definicao de PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// --- Fatores de Conversão ---
#define POLEGADA_PARA_METRO (0.0254)
#define LBF_PARA_NEWTON (4.4482216153)
#define LBIN_PARA_NM (LBF_PARA_NEWTON * POLEGADA_PARA_METRO)
#define PSI_PARA_PASCAL (LBF_PARA_NEWTON / (POLEGADA_PARA_METRO * POLEGADA_PARA_METRO))
#define PASCAL_PARA_MPA (1.0 / 1000000.0)

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
#define Q1_ANGULO_FORCA_POLIA_D_DEG 40.0
#define Q1_FATOR_FORCA_POLIA_V 1.5

// Posicoes (polegadas) - Mancais em A e C
#define Q1_X_A 0.0
#define Q1_X_B 10.0 // Engrenagem B
#define Q1_X_D 20.0 // Polia D
#define Q1_X_C 26.0 // Mancal C

// --- Propriedades do Material (Aco SAE 1040 Estirado a Frio) e Projeto ---
#define Q1_S_U_PSI 85000.0       // Resistencia a tracao ultima (psi) para SAE 1040 CD
#define Q1_S_Y_PSI 71000.0       // Tensao de escoamento (psi) para SAE 1040 CD
#define Q1_S_N_BASE_PSI 38000.0  // Limite de fadiga basico (psi) para SAE 1040 CD, usinada/trefilada (Fig 5.8 [1])
#define Q1_C_SIZE_DEFAULT 0.85   // Fator de tamanho (estimativa inicial)
#define Q1_C_RELIAB_99 0.81      // Fator de confiabilidade para 99%
#define Q1_FATOR_PROJETO_N 2.5   // Fator de seguranca
// Kt para pontos de interesse (assumindo rasgos de chaveta e filetes de mancal)
#define Q1_KT_CHAVETA 2.0
#define Q1_KT_MANCAL_FILETE 2.5


// Estrutura para armazenar forcas em um ponto (Fy, Fz em lb)
typedef struct {
    double Fy_lb;
    double Fz_lb;
} Q1_ForcaImperial;

