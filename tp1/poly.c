#include <stdio.h>
#include <stdlib.h>

#include "poly.h"

#include <x86intrin.h>

p_polyf_t creer_polynome(int degre) {
  p_polyf_t p;

  p = (p_polyf_t)malloc(sizeof(polyf_t));
  p->degre = degre;

  p->coeff = (float *)malloc((degre + 1) * sizeof(float));

  return p;
}

void detruire_polynome(p_polyf_t p) {
  free(p->coeff);
  free(p);

  return;
}

void init_polynome(p_polyf_t p, float x) {
  register unsigned int i;

  for (i = 0; i <= p->degre; ++i)
    p->coeff[i] = x;

  return;
}

p_polyf_t lire_polynome_float(char *nom_fichier) {
  FILE *f;
  p_polyf_t p;
  int degre;
  int i;
  int cr;

  f = fopen(nom_fichier, "r");
  if (f == NULL) {
    fprintf(stderr, "erreur ouverture %s \n", nom_fichier);
    exit(-1);
  }

  cr = fscanf(f, "%d", &degre);
  if (cr != 1) {
    fprintf(stderr, "erreur lecture du degre\n");
    exit(-1);
  }
  p = creer_polynome(degre);

  for (i = 0; i <= degre; i++) {
    cr = fscanf(f, "%f", &p->coeff[i]);
    if (cr != 1) {
      fprintf(stderr, "erreur lecture coefficient %d\n", i);
      exit(-1);
    }
  }

  fclose(f);

  return p;
}

void ecrire_polynome_float(p_polyf_t p) {
  int i;

  printf("%f + %f x ", p->coeff[0], p->coeff[1]);

  for (i = 2; i <= p->degre; i++) {
    printf("+ %f X^%d ", p->coeff[i], i);
  }

  printf("\n");

  return;
}

#define EPSILON 0.00000001

int egalite_polynome(p_polyf_t p1, p_polyf_t p2) {
  if (p1->degre == p2->degre) {
    for (int i = 0; i < p1->degre; i++) {
      if (p1->coeff[i] - p2->coeff[i] < EPSILON) {
        return 0;
      }
    }
    return 1;
  }
  return 0;
}

p_polyf_t addition_polynome(p_polyf_t p1, p_polyf_t p2) {
  p_polyf_t p3;
  register unsigned int i;

  p3 = creer_polynome(max(p1->degre, p2->degre));

  for (i = 0; i <= min(p1->degre, p2->degre); ++i) {
    p3->coeff[i] = p1->coeff[i] + p2->coeff[i];
  }

  if (p1->degre > p2->degre) {
    for (i = (p2->degre + 1); i <= p1->degre; ++i)
      p3->coeff[i] = p1->coeff[i];
  } else if (p2->degre > p1->degre) {
    for (i = (p1->degre + 1); i <= p2->degre; ++i)
      p3->coeff[i] = p2->coeff[i];
  }

  return p3;
}

p_polyf_t multiplication_polynome_scalaire(p_polyf_t p, float alpha) {
  p_polyf_t p2 = creer_polynome(p->degre);
  for (int i = 0; i < p->degre; i++) {
    p2->coeff[i] = alpha * p->coeff[i];
  }
  return p2;
}

float eval_polynome(p_polyf_t p, float x) {
  float res = 0.0;
  float x_puiss = 1.0;
  for (int i = 0; i < p->degre; i++) {
    res += p->coeff[i] * x_puiss;
    x_puiss *= x;
  }
  return res;
}

p_polyf_t multiplication_polynomes(p_polyf_t p1, p_polyf_t p2) {
  if (p1->degre < p2->degre) {
    return multiplication_polynomes(p2, p1);
  } else {
    p_polyf_t res = creer_polynome(p2->degre);
    for (int i = 0; i < p2->degre; i++) {
      p_polyf_t p1_fois_p2i = creer_polynome(p1->degre);
      for (int j = 0; j < p1->degre; j++) {
        p1_fois_p2i->coeff[j] = p2->coeff[i] * p1->coeff[j];
      }
      res = addition_polynome(res, p1_fois_p2i);
    }
    return res;
  }
}

p_polyf_t puissance_polynome(p_polyf_t p, int n) {
  p_polyf_t res = creer_polynome(p->degre + n);
  for (int i = n; i < res->degre; i++) {
    res[i] = p[i - n];
  }

  return res;
}

p_polyf_t composition_polynome(p_polyf_t p, p_polyf_t q) {
  p_polyf_t res = creer_polynome(p->degre * q->degre);
  for (int i = 0; i < p->degre; i++) {
    res = addition_polynome(res, multiplication_polynome_scalaire(
                                     puissance_polynome(q, i), p->coeff[i]));
  }
  return res;
}

p_polyf_creux_t init_poly_creux(int degre) {
    p_polyf_creux_t p = malloc(sizeof(polyf_creux_t));
    p->coeff = 0.0;
    p->degre = degre;
    p->suivant = NULL;
    return p;
}

void detruire_poly_creux(p_polyf_creux_t p) {
    if (p != NULL) {
        detruire_poly_creux(p->suivant);
        free(p);
    }
}

p_polyf_creux_t addition_poly_creux(p_polyf_creux_t p, p_polyf_creux_t q) {
    p_polyf_creux_t res = NULL;
    while (p) {
        res->suivant = init_poly_creux(p->degre);
        p_polyf_creux_t q2 = q;
        res->coeff = p->coeff;
        p = p->suivant;
    }
    while (q) {
        p_polyf_creux_t p2 = res;
        int egal = 0;
        while (p2 && !egal) {
            if (p2->degre == q->degre) {
                p2->coeff += q->coeff;
                egal = 1;
            }
            p2 = p2->suivant;
        }
        q = q->suivant;
    }
    return res;
}

p_polyf_creux_t multiplication_poly_creux_scalaire(p_polyf_creux_t p, float s) {
    p_polyf_creux_t res = NULL;
    while (p) {
        p_polyf_creux_t monome = init_poly_creux(p->degre);
        monome->coeff = s * p->coeff;
        monome->suivant = res;
        res = monome;
    }
    return res;
}
