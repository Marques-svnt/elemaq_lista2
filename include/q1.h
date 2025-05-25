#ifndef Q1_H
#define Q1_H

void q1(FILE *fp_respostas, FILE *fp_dados);

// Definicao de PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// --- Constantes de entrada para a Questao 1 ---
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
#define Q1_FATOR_FORCA_POLIA_V 1.5      // F_total = Fator * F_efetiva

// Posicoes (polegadas) - Mancais em A e C (do problema original Q1)
#define X_A 0.0
#define X_B 10.0 // Engrenagem B
#define X_D 20.0 // Polia D
#define X_C 26.0 // Mancal C (final do eixo considerado para Q1)

// Estrutura para armazenar forcas em um ponto
typedef struct {
    double Fy;
    double Fz;
} Forca;


#endif 