#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "q1.h"
#include "q2.h"
#include "q3.h"
#include "q4.h"

int main()
{
    FILE *arquivo_respostas;
    FILE *arquivo_dados;
    arquivo_respostas = fopen("respostas.txt", "w");
    arquivo_dados = fopen("dados.txt", "w");

    q1(arquivo_respostas, arquivo_dados);

    q2(arquivo_respostas, arquivo_dados);

    q3(arquivo_respostas, arquivo_dados);

    q4(arquivo_respostas, arquivo_dados);

    fclose(arquivo_respostas);
    fclose(arquivo_dados);

    int status_diagram = system("python src/diagram.py");
    if (status_diagram != 0)
    { // Simplificando, qualquer status diferente de 0 é um problema
        fprintf(stderr, "O script Python retornou um erro (status: %d).\n", status_diagram);
        fprintf(stderr, "Isso pode ser devido a dependencias Python ausentes (matplotlib, numpy).\n");
        fprintf(stderr, "Por favor, certifique-se de que Python 3 esta instalado e execute:\n");
        fprintf(stderr, "pip install matplotlib numpy\n");
        fprintf(stderr, "Ou execute o setup.bat\n");
    }

    int status_pdf = system("python src/pdf.py");
    if (status_pdf != 0)
    { // Simplificando, qualquer status diferente de 0 é um problema
        fprintf(stderr, "O script Python retornou um erro (status: %d).\n", status_pdf);
        fprintf(stderr, "Isso pode ser devido a dependencias Python ausentes (fpdf2, pillow).\n");
        fprintf(stderr, "Por favor, certifique-se de que Python 3 esta instalado e execute:\n");
        fprintf(stderr, "pip install fpdf2 pillow\n");
        fprintf(stderr, "Ou use o arquivo: pip install -r requirements.txt\n");
        fprintf(stderr, "Ou execute o setup.bat\n");
    }
    return 0;
}