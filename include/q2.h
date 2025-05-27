#ifndef Q2_H
#define Q2_H

void q2(FILE *fp_respostas, FILE *fp_dados);

// Definicao de PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// --- Fatores de Conversão ---
#define POLEGADA_PARA_METRO (0.0254)
#define LBF_PARA_NEWTON (4.4482216153)
#define LBIN_PARA_NM (LBF_PARA_NEWTON * POLEGADA_PARA_METRO)

// --- Constantes de entrada para a Questao 2 (Imperial) ---
#define Q2_POTENCIA_HP 20.0
#define Q2_ROTACAO_RPM 750.0
#define Q2_FATOR_TORQUE 63000.0 // lb.in * rpm / hp

// Parametros da Engrenagem B (Questao 2)
#define Q2_PASSO_DIAMETRAL_GEARB 6.0
#define Q2_NUM_DENTES_GEARB 100.0
#define Q2_ANGULO_PRESSAO_GEARB_DEG 20.0 // Mantido da Q1

// Parametros da Polia D (Questao 2)
#define Q2_DIAMETRO_POLIA_D_IN 9.0
#define Q2_ANGULO_FORCA_POLIA_D_DEG 40.0 // Angulo com a vertical (Mantido da Q1)
#define Q2_FATOR_FORCA_POLIA_V 1.5

// Posicoes (polegadas) - Mantidas da estrutura da Q1
#define Q2_X_A 0.0  // Mancal A
#define Q2_X_B 10.0 // Engrenagem B
#define Q2_X_D 20.0 // Polia D
#define Q2_X_C 26.0 // Mancal C (final do eixo considerado para Q1/Q2)

// Estrutura para armazenar forcas em um ponto (em lb)
typedef struct {
    double Fy_lb;
    double Fz_lb;
} Q2_ForcaImperial;


#endif 