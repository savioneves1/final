#include <stdio.h>
#include <stdlib.h>

struct tipofila {
    char dado;
    struct tipofila *prox;
} ;

void InicializaFila (struct tipofila **fila) {
    *fila = NULL;
    return;
}

int FilaVazia (struct tipofila *fila) {
    if (fila == NULL)
        return 1;
    else
        return 0;
}

void InsereFila(struct tipofila **fila, char dadonovo) {
    struct tipofila *f1, *f2;
    f1 = malloc (sizeof (struct tipofila));
    f1->dado = dadonovo;
    f1->prox = NULL;
    if (*fila == NULL)
        *fila = f1;
    else {
        f2 = *fila;
        while (f2->prox != NULL)
            f2 = f2->prox;
        f2->prox = f1;
    }
    return;
}

char RetiraFila(struct tipofila **fila) {
    struct tipofila *f1;
    char car;
    f1 = *fila;
    *fila = f1->prox;
    car = f1->dado;
    free (f1);
    return car;
}

int main(){

    FILE *arq, *arqT, *arqC;
    char c;
    struct tipofila *finicio;
    int count=0;

    int count1=42;

    InicializaFila(&finicio);

    arq = fopen ("dadosdctphoto.txt", "r");
    arqT = fopen ("dadosdctphoto_temp.txt", "w");
    arqC = fopen ("dadosdctphoto_coef.txt", "w");

    while ((c = getc (arq)) != EOF) {
        if (c!='\n'){
            count=count+1;
            InsereFila(&finicio,c);
        }
        if(c=='\n'){
            InsereFila(&finicio,c);
            if(count==count1){
                while(!FilaVazia(finicio)){
                    fputc(RetiraFila(&finicio),arqC);
                }
            }
            else{
                while(!FilaVazia(finicio)){
                    fputc(RetiraFila(&finicio),arqT);
                }
            }
            count=0;
        }
    }
    fclose (arq);
    fclose (arqT);
    fclose(arqC);
    return 0;
}
