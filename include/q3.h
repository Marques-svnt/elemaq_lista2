#ifndef Q3_H
#define Q3_H

void q3(FILE *fp_respostas, FILE *fp_dados);

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// --- Constantes da Questao 3 ---
#define Q3_VELOCIDADE_RPM 200.0
#define Q3_FATOR_TORQUE 63000.0

// Polia A (Correia Plana)
#define Q3_POTENCIA_POLIA_A_HP 10.0
#define Q3_DIAMETRO_POLIA_A_IN 20.0
#define Q3_RELACAO_F1_F2_POLIA_A 2.5

// Engrenagem C (Dentes Retos)
#define Q3_POTENCIA_ENGRENAGEM_C_HP 6.0
#define Q3_PASSO_DIAMETRAL_ENG_C 8.0
#define Q3_NUM_DENTES_ENG_C 80.0
#define Q3_ANGULO_PRESSAO_ENG_C_DEG 20.0

// Roda Dentada D (Corrente)
#define Q3_POTENCIA_RODA_DENTADA_D_HP 4.0
#define Q3_DIAMETRO_RODA_DENTADA_D_IN 6.0

// Posicoes dos Elementos e Mancais (polegadas)
#define Q3_X_A 0.0  // Polia A (Entrada)
#define Q3_X_B 5.0  // Mancal B
#define Q3_X_C 13.0 // Engrenagem C (Saida)
#define Q3_X_D 21.0 // Roda Dentada D (Saida)
#define Q3_X_E 26.0 // Mancal E

// Propriedades do Material (para respostas.txt, não para dados.txt de diagramas)
#define Q3_S_U_PSI 69000.0
#define Q3_S_Y_PSI 58000.0
#define Q3_S_N_SURF_PSI 28000.0
#define Q3_C_RELIAB_99 0.81
#define Q3_FATOR_PROJETO_N 2.5
#define Q3_KT_POLIA_A 2.0
#define Q3_KT_MANCAL_B 2.5
#define Q3_KT_ENGRENAGEM_C 2.0
#define Q3_KT_RODA_DENTADA_D 2.0
#define Q3_KT_MANCAL_E 2.5

// Estrutura para forcas (reutilizada)
typedef struct {
    double Fy;
    double Fz;
} Q3_ForcaPonto;



#endif 