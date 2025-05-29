#ifndef Q2_H
#define Q2_H

void q2(FILE *fp_respostas, FILE *fp_dados);

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// --- Fatores de Conversão ---
#define POLEGADA_PARA_METRO (0.0254)
#define LBF_PARA_NEWTON (4.4482216153)
#define LBIN_PARA_NM (LBF_PARA_NEWTON * POLEGADA_PARA_METRO)
#define PSI_PARA_PASCAL (LBF_PARA_NEWTON / (POLEGADA_PARA_METRO * POLEGADA_PARA_METRO))
#define PASCAL_PARA_MPA (1.0 / 1000000.0)

// --- Constantes de entrada para a Questao 2 (Imperial) ---
#define Q2_POTENCIA_HP 20.0
#define Q2_ROTACAO_RPM 750.0
#define Q2_FATOR_TORQUE 63000.0

// Engrenagem B (Questao 2)
#define Q2_PASSO_DIAMETRAL_GEARB 6.0
#define Q2_NUM_DENTES_GEARB 100.0
#define Q2_ANGULO_PRESSAO_GEARB_DEG 20.0

// Polia D (Questao 2)
#define Q2_DIAMETRO_POLIA_D_IN 9.0
#define Q2_ANGULO_FORCA_POLIA_D_DEG 40.0
#define Q2_FATOR_FORCA_POLIA_V 1.5

// Posicoes (polegadas) - Mantidas da estrutura da Q1
#define Q2_X_A 0.0
#define Q2_X_B 10.0
#define Q2_X_D 20.0
#define Q2_X_C 26.0

// --- Propriedades do Material (Aço SAE 1040 Estirado a Frio - Imperial) ---
#define Q2_S_U_PSI 85000.0      // Resistencia a Tracao Ultima (psi)
#define Q2_S_Y_PSI 71000.0      // Tensao de Escoamento (psi)
#define Q2_S_N_BASE_PSI 39000.0 // Limite de Fadiga Base Sn (para usinada/estirada a frio, de Su=85ksi e Fig. 5.8) [cite: 1]

// --- Parametros de Projeto para Cálculo de Diametro ---
#define Q2_C_SIZE_PADRAO 0.85   // Fator de Tamanho Cs inicial (Fig. 5.9) [cite: 1]
#define Q2_C_RELIAB_99 0.81     // Fator de Confiabilidade CR para 99% (Tabela 5.2) [cite: 1]
#define Q2_FATOR_PROJETO_N 2.5  // Fator de Seguranca N

// Fatores de Concentracao de Tensao (Kt) Adimensionais - Assumidos
#define Q2_KT_MANCAL 2.5        // Para filetes agudos/degraus em mancais
#define Q2_KT_ENG_POLIA 2.2     // Para rasgos de chaveta em engrenagem/polia

// Estrutura para armazenar forcas em um ponto (em lb)
typedef struct {
    double Fy_lb;
    double Fz_lb;
} Q2_ForcaImperial;



#endif 