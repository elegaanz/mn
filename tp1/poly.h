/*
  polyf_t   : structure polynome
  p_polyf_t : pointeur sur un polynome
*/

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

typedef struct {
  int degre;
  float *coeff;
} polyf_t, *p_polyf_t;

p_polyf_t creer_polynome(int degre);

void init_polynome(p_polyf_t p, float x);

void detruire_polynome(p_polyf_t p);

p_polyf_t lire_polynome_float(char *nom_fichier);

void ecrire_polynome_float(p_polyf_t p);

int egalite_polynome(p_polyf_t p1, p_polyf_t p2);

p_polyf_t addition_polynome(p_polyf_t p1, p_polyf_t p2);

p_polyf_t multiplication_polynome_scalaire(p_polyf_t p, float alpha);

float eval_polynome(p_polyf_t p, float x);

p_polyf_t multiplication_polynomes(p_polyf_t p1, p_polyf_t p2);

p_polyf_t puissance_polynome(p_polyf_t p, int n);

p_polyf_t composition_polynome(p_polyf_t p, p_polyf_t q);

typedef struct poly_creux {
  int degre;
  float coeff;
  poly_creux *suivant;
} polyf_creux_t, *p_polyf_creux_t;

p_polyf_creux_t init_poly_creux(int degre);

void detruire_poly_creux(p_polyf_creux_t p);

p_polyf_creux_t addition_poly_creux(p_polyf_creux_t p, p_polyf_creux_t q);

p_polyf_creux_t multiplication_poly_creux_scalaire(p_polyf_creux_t, float s);
