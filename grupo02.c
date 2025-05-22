
#include <stdio.h>
#include <stdlib.h> // Para fabs
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Estrutura para um ponto em um diagrama de esforços
typedef struct
{
    double x;     // Posição ao longo do eixo (m)
    double valor; // Valor do esforço (N para cortante, N-m para momento)
} PontoDiagrama;

// Estrutura para dados de um ponto crítico para cálculo de diâmetro
typedef struct
{
    double x_pos;                     // Posição do ponto crítico (m)
    double momento_fletor_resultante; // Momento fletor M (N-m)
    double momento_torsor;            // Momento torsor T (N-m)
} PontoCriticoParaDiametro;

// --- Módulo 2.1: Geração de Dados para Diagramas de Esforços ---
int gerar_pontos_diagrama_cortante(
    double comprimento_eixo,
    int num_cargas_verticais,
    const double pos_cargas_verticais[],
    const double val_cargas_verticais[],
    PontoDiagrama pontos_cortante[],
    int max_pontos_cortante)
{

    int contador_pontos = 0;
    if (contador_pontos < max_pontos_cortante)
    {
        pontos_cortante[contador_pontos++] = (PontoDiagrama){0.0, 0.0};
    }

    if (num_cargas_verticais > 0 && contador_pontos < max_pontos_cortante - 1)
    {
        double cortante_atual = pontos_cortante[contador_pontos - 1].valor; // Cortante antes da carga
        if (pos_cargas_verticais[0] > pontos_cortante[contador_pontos - 1].x)
        { // Se a carga não está em x=0
            if (contador_pontos < max_pontos_cortante)
            {
                pontos_cortante[contador_pontos++] = (PontoDiagrama){pos_cargas_verticais[0], cortante_atual};
            }
        }
        cortante_atual += val_cargas_verticais[0]; // Cortante após a carga
        if (contador_pontos < max_pontos_cortante)
        {
            // Se a posição da carga já foi adicionada, atualiza o valor do último ponto ou adiciona novo
            if (contador_pontos > 0 && pontos_cortante[contador_pontos - 1].x == pos_cargas_verticais[0])
            {
                pontos_cortante[contador_pontos - 1].valor = cortante_atual; // Atualiza se x é o mesmo
            }
            else
            {
                pontos_cortante[contador_pontos++] = (PontoDiagrama){pos_cargas_verticais[0], cortante_atual};
            }
        }
    }
    // *** FIM DA LÓGICA DE ESTÁTICA ***

    // Exemplo final: Adiciona ponto no final do eixo, idealmente fechando em 0
    if (contador_pontos > 0 && contador_pontos < max_pontos_cortante)
    {
        double ultimo_cortante = pontos_cortante[contador_pontos - 1].valor;
        if (pontos_cortante[contador_pontos - 1].x < comprimento_eixo)
        {
            pontos_cortante[contador_pontos++] = (PontoDiagrama){comprimento_eixo, ultimo_cortante};
        }
        // Verifica se o cortante no final do eixo é zero (para sistemas em equilíbrio)
        if (fabs(ultimo_cortante) > 1e-6 && contador_pontos < max_pontos_cortante)
        {
            pontos_cortante[contador_pontos++] = (PontoDiagrama){comprimento_eixo, 0.0};
        }
    }
    return contador_pontos;
}

int gerar_pontos_diagrama_fletor(
    double comprimento_eixo,
    int num_pontos_cortante,
    const PontoDiagrama pontos_cortante[], // Deve estar ordenado por x
    int num_momentos_concentrados,
    const double pos_momentos_concentrados[], // Deve estar ordenado por x
    const double val_momentos_concentrados[],
    PontoDiagrama pontos_fletor[],
    int max_pontos_fletor)
{

    int contador_pontos = 0;
    double momento_atual = 0.0;
    int idx_mom_conc = 0; // Índice para percorrer momentos concentrados

    if (contador_pontos < max_pontos_fletor)
    {
        pontos_fletor[contador_pontos++] = (PontoDiagrama){0.0, 0.0};
    }

    // *** INÍCIO DA LÓGICA DE INTEGRAÇÃO E MOMENTOS (A SER IMPLEMENTADA PELO GRUPO) ***
    for (int i = 0; i < num_pontos_cortante - 1 && contador_pontos < max_pontos_fletor; ++i)
    {
        // Considera o intervalo entre pontos_cortante[i] e pontos_cortante[i+1]
        double x_inicio_intervalo = pontos_cortante[i].x;
        double x_fim_intervalo = pontos_cortante[i + 1].x;
        double cortante_no_intervalo = pontos_cortante[i].valor; // Assume cortante constante no intervalo

        if (x_fim_intervalo <= x_inicio_intervalo)
            continue; // Pula se não há avanço em x

        double proximo_x_interesse = x_fim_intervalo;

        // Verifica se há momento concentrado dentro ou no final deste intervalo de cortante
        if (idx_mom_conc < num_momentos_concentrados && pos_momentos_concentrados[idx_mom_conc] < x_fim_intervalo)
        {
            proximo_x_interesse = pos_momentos_concentrados[idx_mom_conc];
        }

        // Calcula momento até proximo_x_interesse
        momento_atual += cortante_no_intervalo * (proximo_x_interesse - x_inicio_intervalo);
        if (contador_pontos < max_pontos_fletor)
        {
            pontos_fletor[contador_pontos++] = (PontoDiagrama){proximo_x_interesse, momento_atual};
        }

        // Aplica momento concentrado, se houver, em proximo_x_interesse
        if (idx_mom_conc < num_momentos_concentrados && fabs(proximo_x_interesse - pos_momentos_concentrados[idx_mom_conc]) < 1e-6)
        {
            momento_atual += val_momentos_concentrados[idx_mom_conc];
            // Atualiza o valor do momento no mesmo ponto x, se já adicionado
            if (contador_pontos > 0 && pontos_fletor[contador_pontos - 1].x == proximo_x_interesse)
            {
                pontos_fletor[contador_pontos - 1].valor = momento_atual;
            }
            else if (contador_pontos < max_pontos_fletor)
            { // Adiciona novo ponto para o salto
                pontos_fletor[contador_pontos++] = (PontoDiagrama){proximo_x_interesse, momento_atual};
            }
            idx_mom_conc++;
        }

        // Se o proximo_x_interesse era um momento concentrado antes do fim do intervalo de cortante constante
        if (proximo_x_interesse < x_fim_intervalo && contador_pontos < max_pontos_fletor)
        {
            momento_atual += cortante_no_intervalo * (x_fim_intervalo - proximo_x_interesse);
            pontos_fletor[contador_pontos++] = (PontoDiagrama){x_fim_intervalo, momento_atual};
        }
    }
    // *** FIM DA LÓGICA DE INTEGRAÇÃO E MOMENTOS ***

    // Garante que o último ponto seja o final do eixo
    if (contador_pontos > 0 && pontos_fletor[contador_pontos - 1].x < comprimento_eixo && contador_pontos < max_pontos_fletor)
    {
        pontos_fletor[contador_pontos++] = (PontoDiagrama){comprimento_eixo, momento_atual}; // ou 0.0
    }

    return contador_pontos;
}
int gerar_pontos_diagrama_torsor(
    double comprimento_eixo,
    int num_torques_aplicados,
    const double pos_torques[], // Deve estar ordenado por x
    const double val_torques[],
    PontoDiagrama pontos_torsor[],
    int max_pontos_torsor)
{

    int contador_pontos = 0;
    double torsor_atual = 0.0;

    if (contador_pontos < max_pontos_torsor)
    {
        pontos_torsor[contador_pontos++] = (PontoDiagrama){0.0, 0.0};
    }

    // *** INÍCIO DA LÓGICA DE TORQUES (A SER IMPLEMENTADA PELO GRUPO) ***
    for (int i = 0; i < num_torques_aplicados && contador_pontos < max_pontos_torsor - 1; ++i)
    {
        // Ponto antes da aplicação do torque i, se a posição for diferente da anterior
        if (pos_torques[i] > pontos_torsor[contador_pontos - 1].x)
        {
            pontos_torsor[contador_pontos++] = (PontoDiagrama){pos_torques[i], torsor_atual};
        }
        // Ponto após a aplicação do torque i
        torsor_atual += val_torques[i];
        // Se a posição da carga já foi adicionada, atualiza o valor do último ponto ou adiciona novo
        if (contador_pontos > 0 && pontos_torsor[contador_pontos - 1].x == pos_torques[i])
        {
            pontos_torsor[contador_pontos - 1].valor = torsor_atual;
        }
        else if (contador_pontos < max_pontos_torsor)
        {
            pontos_torsor[contador_pontos++] = (PontoDiagrama){pos_torques[i], torsor_atual};
        }
    }
    // *** FIM DA LÓGICA DE TORQUES ***

    if (contador_pontos > 0 && pontos_torsor[contador_pontos - 1].x < comprimento_eixo && contador_pontos < max_pontos_torsor)
    {
        pontos_torsor[contador_pontos++] = (PontoDiagrama){comprimento_eixo, torsor_atual};
    }
    // Verifica se o torque no final do eixo é zero (para sistemas em equilíbrio torsional)
    if (contador_pontos > 0 && fabs(pontos_torsor[contador_pontos - 1].x - comprimento_eixo) < 1e-6)
    {
        if (fabs(torsor_atual) > 1e-6 && contador_pontos < max_pontos_torsor)
        {
            pontos_torsor[contador_pontos++] = (PontoDiagrama){comprimento_eixo, 0.0};
        }
    }
    return contador_pontos;
}

// --- Módulo 2.2: Análise de Tensão e Dimensionamento do Eixo por Resistência ---

double calcular_diametro_eixo_fadiga(
    double momento_fletor_M,    // N-m
    double momento_torsor_T,    // N-m
    double Kt_flexao,           // Fator de concentração de tensão para flexão
    double Sy_escoamento,       // Tensão de escoamento (Pa)
    double Sn_fadiga_corrigida, // Resistência à fadiga corrigida (Pa)
    double N_fator_seguranca)
{

    if (Sn_fadiga_corrigida <= 1e-9 || Sy_escoamento <= 1e-9 || N_fator_seguranca <= 1e-9)
    {
        fprintf(stderr, "ERRO: Parâmetros inválidos para cálculo de diâmetro (Sy, Sn', N devem ser > 0).\n");
        return -1.0;
    }

    momento_fletor_M = fabs(momento_fletor_M);
    momento_torsor_T = fabs(momento_torsor_T);

    double termo_flexao_quad = pow((Kt_flexao * momento_fletor_M) / Sn_fadiga_corrigida, 2.0);
    double termo_torsao_quad = (3.0 / 4.0) * pow(momento_torsor_T / Sy_escoamento, 2.0);

    double dentro_raiz = termo_flexao_quad + termo_torsao_quad;
    if (dentro_raiz < 0)
    {
        fprintf(stderr, "ERRO: Valor negativo dentro da raiz no cálculo do diâmetro.\n");
        return -1.0;
    }

    double diametro_cubado = (32.0 * N_fator_seguranca / M_PI) * sqrt(dentro_raiz);
    if (diametro_cubado < 0)
    {
        fprintf(stderr, "ERRO: Diâmetro ao cubo resultou em valor negativo.\n");
        return -1.0;
    }
    return pow(diametro_cubado, 1.0 / 3.0); // Retorna diâmetro em metros
}

// --- Função Principal e Saída de Dados ---
void escrever_dados_para_python(
    const char *nome_arquivo,
    int num_cv, const PontoDiagrama Pcv[],
    int num_fv, const PontoDiagrama Pfv[],
    int num_t, const PontoDiagrama Pt[],
    double comprimento_eixo)
{

    FILE *arquivo = fopen(nome_arquivo, "w");
    if (arquivo == NULL)
    {
        perror("ERRO ao abrir arquivo para escrita");
        return;
    }

    fprintf(arquivo, "# Dados dos Diagramas de Esforços do Eixo (Comprimento: %.3f m)\n", comprimento_eixo);
    fprintf(arquivo, "# Formato: Posicao(m) Torsor(N-m) Fletor(N-m) Cortante(N)\n");

    double val_cv = 0.0, val_fv = 0.0, val_t = 0.0;
    int num_passos_plot = 200; // Número de pontos para amostrar para o arquivo de plotagem

    for (int i = 0; i <= num_passos_plot; ++i)
    {
        double x_atual = (comprimento_eixo / num_passos_plot) * i;

        // Busca valor de cortante
        // Reinicia a busca do índice para cada x_atual para garantir corretude na interpolação degrau
        val_cv = 0.0; // Valor padrão se nenhum ponto encontrado antes de x_atual
        if (num_cv > 0)
            val_cv = Pcv[0].valor; // Assume valor inicial se x_atual for menor que o primeiro ponto
        for (int k = 0; k < num_cv; ++k)
        {
            if (Pcv[k].x <= x_atual)
            {
                val_cv = Pcv[k].valor;
            }
            else
            {
                break;
            }
        }

        // Busca valor de fletor
        val_fv = 0.0;
        if (num_fv > 0)
            val_fv = Pfv[0].valor;
        for (int k = 0; k < num_fv; ++k)
        {
            if (Pfv[k].x <= x_atual)
            {
                val_fv = Pfv[k].valor;
            }
            else
            {
                break;
            }
        }

        // Busca valor de torsor
        val_t = 0.0;
        if (num_t > 0)
            val_t = Pt[0].valor;
        for (int k = 0; k < num_t; ++k)
        {
            if (Pt[k].x <= x_atual)
            {
                val_t = Pt[k].valor;
            }
            else
            {
                break;
            }
        }
        fprintf(arquivo, "%.4f %.2f %.2f %.2f\n", x_atual, val_t, val_fv, val_cv);
    }

    fclose(arquivo);
    printf("Dados dos diagramas escritos em %s (com amostragem simplificada).\n", nome_arquivo);
}

int main()
{
    printf("Simulador de Análise de Eixo - Grupo 2\n");

    // --- Entradas (Exemplo - devem vir de outro módulo ou problema específico) ---
    double comprimento_total_eixo = 1.0; // m

    // Dados para Diagrama de Cortante e Fletor (Exemplo)
    // Reações e Cargas: {posição(m), valor(N)}
    // Ex: Viga biapoiada com carga no meio. Apoios em 0 e 1m. Carga de -2000N em 0.5m.
    // Reações R_A = 1000N, R_B = 1000N.
    double pos_cargas_v[] = {0.0, 0.5, 1.0};
    double val_cargas_v[] = {1000.0, -2000.0, 1000.0}; // R_A, Carga P, R_B
    int num_cargas_v = sizeof(pos_cargas_v) / sizeof(pos_cargas_v[0]);

    // Momentos Concentrados: {posição(m), valor(N-m)}
    double pos_momentos_conc_v[] = {}; // Sem momentos concentrados neste exemplo
    double val_momentos_conc_v[] = {};
    int num_momentos_conc_v = 0;

    // Dados para Diagrama de Momento Torsor (Exemplo)
    // Torques Aplicados: {posição(m), valor(N-m)}
    double pos_torques_aplicados[] = {0.2, 0.8};
    double val_torques_aplicados[] = {150.0, -150.0};
    int num_torques_aplicados = sizeof(pos_torques_aplicados) / sizeof(pos_torques_aplicados[0]);

#define MAX_PONTOS_DIAGRAMA 100 // Aumentado para acomodar mais pontos da lógica de estática
    PontoDiagrama diagrama_cortante[MAX_PONTOS_DIAGRAMA];
    PontoDiagrama diagrama_fletor[MAX_PONTOS_DIAGRAMA];
    PontoDiagrama diagrama_torsor[MAX_PONTOS_DIAGRAMA];

    int n_pontos_cv = 0;
    int n_pontos_fv = 0;
    int n_pontos_t = 0;

    // --- Módulo 2.1: Geração de Dados para Diagramas ---
    // As implementações destas funções são cruciais.
    n_pontos_cv = gerar_pontos_diagrama_cortante(comprimento_total_eixo, num_cargas_v, pos_cargas_v,
                                                 val_cargas_v, // Argumento adicionado
                                                 diagrama_cortante, MAX_PONTOS_DIAGRAMA);
    n_pontos_fv = gerar_pontos_diagrama_fletor(comprimento_total_eixo, n_pontos_cv, diagrama_cortante,
                                               num_momentos_conc_v, pos_momentos_conc_v, val_momentos_conc_v,
                                               diagrama_fletor, MAX_PONTOS_DIAGRAMA);
    n_pontos_t = gerar_pontos_diagrama_torsor(comprimento_total_eixo, num_torques_aplicados, pos_torques_aplicados, val_torques_aplicados,
                                              diagrama_torsor, MAX_PONTOS_DIAGRAMA);

    escrever_dados_para_python("dados_diagramas.txt",
                               n_pontos_cv, diagrama_cortante,
                               n_pontos_fv, diagrama_fletor,
                               n_pontos_t, diagrama_torsor,
                               comprimento_total_eixo);

    // --- Módulo 2.2: Dimensionamento do Eixo por Resistência ---
    printf("\nAnálise de Tensão e Dimensionamento do Eixo:\n");

    double Sy_material = 350e6; // Tensão de escoamento (Pa)
    double Su_material = 550e6; // Tensão de ruptura (Pa)
    double Sn_base_material = 0.5 * Su_material;
    double C_tam = 0.85, C_sup = 0.75, C_conf = 0.814; // Fatores de correção (exemplos)
    double Sn_corrigida_material = Sn_base_material * C_tam * C_sup * C_conf;
    double N_seguranca = 2.0;
    double Kt_f = 1.8; // Fator de concentração de tensão para flexão (ex: rasgo chaveta)

    // Pontos críticos para análise de diâmetro (exemplos)
    // Idealmente, estes seriam determinados a partir da análise dos diagramas gerados
    PontoCriticoParaDiametro pontos_para_analise[] = {
        // {x_pos, MomentoFletorMax, MomentoTorsorNoPonto}
        {0.2, 100.0, 150.0}, // Exemplo: M=100 Nm, T=150 Nm em x=0.2m
        {0.5, 250.0, 150.0}, // Exemplo: M=250 Nm (pico do fletor), T=150 Nm em x=0.5m
        {0.8, 80.0, 0.0}     // Exemplo: M=80 Nm, T=0 Nm (após saída de torque) em x=0.8m
    };
    int num_pontos_para_analise = sizeof(pontos_para_analise) / sizeof(pos_torques_aplicados[0]); // CORREÇÃO: Usar o array correto para sizeof

    printf("Calculando diâmetros mínimos nos pontos críticos (usando Eq. Fadiga):\n");
    for (int i = 0; i < num_pontos_para_analise; ++i)
    {
        double M_analise = pontos_para_analise[i].momento_fletor_resultante;
        double T_analise = pontos_para_analise[i].momento_torsor;

        // Obter M e T dos diagramas gerados para a posição x_pos[i] seria mais preciso.
        // A estrutura PontoCriticoParaDiametro já assume que M e T são conhecidos para aquele ponto.

        double d_minimo = calcular_diametro_eixo_fadiga(M_analise, T_analise, Kt_f,
                                                        Sy_material, Sn_corrigida_material, N_seguranca);

        if (d_minimo > 0)
        {
            printf("Ponto x=%.2f m: M=%.1f Nm, T=%.1f Nm -> D_min = %.4f m (%.2f mm)\n",
                   pontos_para_analise[i].x_pos, M_analise, T_analise, d_minimo, d_minimo * 1000.0);
        }
        else
        {
            printf("Ponto x=%.2f m: Erro no cálculo do diâmetro.\n", pontos_para_analise[i].x_pos);
        }
    }
    return 0;
}
