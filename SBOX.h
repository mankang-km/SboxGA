#ifndef _SBOX_
#define _SBOX_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define POPSIZE 100
#define LENGTH 256
#define BITS 8
#define TOURSIZE 3

#define Pc 0.8
#define Pm 0.3

extern int generation;
extern double sumfitness;

struct individual
{
	int chorm[LENGTH];
	int fitness;
	int duniformity;
	int linearity;
	int buniformity;
};

extern int best;
extern int bestd;
extern int bestb;
extern int bestl;

extern struct individual population[POPSIZE];
extern struct individual newpopulation[POPSIZE];
extern struct individual pool[POPSIZE];
extern struct individual bestp;

extern int qflag;
extern int number;

extern int C[720][6];


extern FILE *fpdata;
extern FILE *fpout;


void GenerateInitPopulation(void);

void CalculateFitnessValue(void);

void SelectionOperation(void);

void CrossoverOperation(void);

void MutationOperation(void);

void cpopulation(int c);

void FindBest(void);



#endif