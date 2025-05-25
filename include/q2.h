#ifndef Q2_H
#define Q2_H

void q2(FILE *fp_respostas, FILE *fp_dados);

// Definicao de PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// --- Constantes de entrada para a Questao 2 ---
#define Q2_POTENCIA_HP 20.0
#define Q2_ROTACAO_RPM 750.0
#define Q2_FATOR_TORQUE 63000.0

// Engrenagem B (Questao 2)
#define Q2_PASSO_DIAMETRAL_GEARB 6.0
#define Q2_NUM_DENTES_GEARB 100.0
#define Q2_ANGULO_PRESSAO_GEARB_DEG 20.0 // Mantido da Q1

// Polia D (Questao 2)
#define Q2_DIAMETRO_POLIA_D_IN 9.0
#define Q2_ANGULO_FORCA_POLIA_D_DEG 40.0 // Angulo com a vertical (Mantido da Q1)
#define Q2_FATOR_FORCA_POLIA_V 1.5      // F_total = Fator * F_efetiva

// Posicoes (polegadas) - Mantidas da estrutura da Q1
#define Q2_X_A 0.0  // Mancal A
#define Q2_X_B 10.0 // Engrenagem B
#define Q2_X_D 20.0 // Polia D
#define Q2_X_C 26.0 // Mancal C (final do eixo considerado para Q1/Q2)

// Estrutura para armazenar forcas em um ponto (reutilizada)
typedef struct {
    double Fy;
    double Fz;
} Q2_Forca;


#endif 