#include "SBOX.h"


void GenerateInitPopulation(void)
{
	int i, j,index,temp;

	for (i = 0; i<POPSIZE; i++)
	{
		for (j = 0; j<LENGTH; j++)
		{
			population[i].chorm[j] = j;
		}
	}

	for (i = 0; i < POPSIZE; i++)
	{
		for (j = 0; j<LENGTH; j++)
		{
			index = rand() % (LENGTH - j);
			temp = population[i].chorm[j];
			population[i].chorm[j] = population[i].chorm[j+ index];
			population[i].chorm[j + index] = temp;
		}

	}

}


int innerproduct(int x, int a)
{
	int sum = 0, index;
	for (index = 0; index < BITS; index++)
	{
		if ((a >> index)&(0x1))
		{
			sum = sum ^ ((x >> index)&(0x1));
		}
	}
	return sum;
}



void CalculateFitnessValue(void)
{
	
	int i, j,k,m,n;
	int DDT[LENGTH][LENGTH], LAT[LENGTH][LENGTH], BCT[LENGTH][LENGTH];
	int detain, detaout, Sin, alpha, beta;
	int sum1, sum2;
	int inverseSbox[LENGTH];
	double p;
	int position, back[4];
	struct individual minf;


	for (i = 0; i < POPSIZE; i++)
	{
		for (j = 0; j < LENGTH; j++)
		{
			for (k = 0; k < LENGTH; k++)
			{
				DDT[j][k] = 0;
				BCT[j][k] = 0;
				LAT[j][k] = -pow((double)2, BITS - 1);
			}
		}

		for (detain = 0; detain < LENGTH; detain++)
		{
			for (Sin = 0; Sin < LENGTH; Sin++)
			{
				detaout = population[i].chorm[Sin] ^ population[i].chorm[Sin^detain];
				DDT[detain][detaout] ++;
			}
		}

		for (alpha = 0; alpha < LENGTH; alpha++)
		{
			for (beta = 0; beta < LENGTH; beta++)
			{
				for (Sin = 0; Sin < LENGTH; Sin++)
				{
					if (innerproduct(Sin, alpha) == innerproduct(population[i].chorm[Sin], beta))
						LAT[alpha][beta] ++;
				}
			}
		}

		for (j = 0; j<LENGTH; j++)
		{
			inverseSbox[population[i].chorm[j]] = j;
		}
		for (detain = 0; detain < LENGTH; detain++)
		{
			for (detaout = 0; detaout< LENGTH; detaout++)
			{
				for (Sin = 0; Sin < LENGTH; Sin++)
				{
					if (detain == ((inverseSbox[population[i].chorm[Sin] ^ detaout]) ^ (inverseSbox[population[i].chorm[Sin^detain] ^ detaout])))
					{
						BCT[detain][detaout] ++;
					}
				}
			}
		}

		population[i].duniformity = 0;
		for (j = 0; j < LENGTH; j++)
		{
			for (k = 0; k < LENGTH; k++)
			{
				if ((DDT[j][k] > population[i].duniformity) && ((j + k) != 0))
					population[i].duniformity = DDT[j][k];
			}
		}
		population[i].linearity = 0;
		for (j = 0; j < LENGTH; j++)
		{
			for (k = 0; k < LENGTH; k++)
			{
				if ((abs(LAT[j][k]) > abs(population[i].linearity)) && ((j + k) != 0))
					population[i].linearity = LAT[j][k];
			}
		}
		population[i].linearity = 2 * population[i].linearity;
		population[i].buniformity = 0;
		for ( j = 0; j < LENGTH; j++)
		{
			for ( k = 0; k < LENGTH; k++)
			{
				if ((BCT[j][k] > population[i].buniformity) && (j != 0) && (k != 0))
					population[i].buniformity = BCT[j][k];
			}
		}

		population[i].fitness = population[i].duniformity + 3*abs(population[i].linearity) + population[i].buniformity;

		
	}

	
}

int comp(const void*a, const void*b)
{
	return (*(struct individual *)a).fitness - (*(struct individual *)b).fitness;
}


void cpopulation(int c)
{
	int i, j;
	qsort(population, POPSIZE, sizeof(population[0]), comp);
	for (i = 0; i < 10; i++)
	{
		pool[c * 10 + i] = population[i];
	}
}

void FindBest(void)
{
	int i, j;

	if (!qflag)
	{
		fpout = fopen("resultout.txt", "w");
		qflag = 1;
	}
	else
	{
		fpout = fopen("resultout.txt", "a");
	}

	for (i = 0; i < POPSIZE; i++)
	{
		if (population[i].fitness < best)
		{
			best = population[i].fitness;
			fprintf(fpout, " the number of iterations is %d \n \n ", generation);
			fprintf(fpout, "find a best SBOX: ");
			for (j = 0; j < LENGTH; j++)
			{
				fprintf(fpout, "%d, ", population[i].chorm[j]);
			}
			fprintf(fpout, "\n The differential uniformity is: %d, The linearity is: %d, The boomerang uniformity is %d \n \n", population[i].duniformity, abs(population[i].linearity), population[i].buniformity);

		}
		if (population[i].duniformity < bestd)
		{
			bestd = population[i].duniformity;
			fprintf(fpout, " the number of iterations is %d \n \n ", generation);
			fprintf(fpout, "find a best differential uniformity SBOX: ");
			for (j = 0; j < LENGTH; j++)
			{
				fprintf(fpout, "%d, ", population[i].chorm[j]);
			}
			fprintf(fpout, "\n The differential uniformity is: %d, The linearity is: %d, The boomerang uniformity is %d \n \n", population[i].duniformity, abs(population[i].linearity), population[i].buniformity);

		}
		if (abs(population[i].linearity) < bestl)
		{
			bestl = abs(population[i].linearity);
			fprintf(fpout, " the number of iterations is %d \n \n ", generation);
			fprintf(fpout, "find a best linearity SBOX: ");
			for (j = 0; j < LENGTH; j++)
			{
				fprintf(fpout, "%d, ", population[i].chorm[j]);
			}
			fprintf(fpout, "\n The differential uniformity is: %d, The linearity is: %d, The boomerang uniformity is %d \n \n", population[i].duniformity, abs(population[i].linearity), population[i].buniformity);

		}
		if (population[i].buniformity < bestb)
		{
			bestb = population[i].buniformity;
			fprintf(fpout, " the number of iterations is %d \n \n ", generation);
			fprintf(fpout, "find a best boomerang uniformity SBOX: ");
			for (j = 0; j < LENGTH; j++)
			{
				fprintf(fpout, "%d, ", population[i].chorm[j]);
			}
			fprintf(fpout, "\n The differential uniformity is: %d, The linearity is: %d, The boomerang uniformity is %d \n \n", population[i].duniformity, abs(population[i].linearity), population[i].buniformity);
		}

	}

	fclose(fpout);
}


void SelectionOperation(void)
{
	int i, j;
	int index1,index2;

	for (i = 0; i < POPSIZE; i++)
	{
		index1 = rand() % POPSIZE;
		for (j = 0; j < TOURSIZE - 1; j++)
		{
			index2 = rand() % POPSIZE;
			if (population[index1].fitness>population[index2].fitness)
			{
				index1 = index2;
			}
		}
		newpopulation[i]= population[index1];
	}

	for (i = 0; i<POPSIZE; i++)
	{
		population[i] = newpopulation[i];
	}
}


void CrossoverOperation(void)
{
	int i, j, k, parent,point[2], temp;
	int index[POPSIZE];
	double p;

	for (i = 0; i<POPSIZE; i++)
	{
		index[i] = i;
	}
	for (i = 0; i<POPSIZE; i++)
	{
		parent = rand() % (POPSIZE - i);
		temp = index[i];
		index[i] = index[i + parent];
		index[i + parent] = temp;
	}

	for (i = 0; i < (POPSIZE - 1); i += 2)
	{
		p = rand() / (double)(RAND_MAX);
		if (p < Pc)
		{
			point[0] = rand() % LENGTH;
			point[1] = rand() % LENGTH;

			if (point[0]>point[1])
			{
				temp = point[0];
				point[0] = point[1];
				point[1] = temp;
			}

			for (j = point[0]; j<point[1]; j++)//½øÐÐ½»²æ
			{
				temp = population[index[i]].chorm[j];
				population[index[i]].chorm[j] = population[index[i + 1]].chorm[j];
				population[index[i + 1]].chorm[j] = temp;
			}

			for (j = 0; j < point[0];j++)
			{
				for (k = point[0]; k < point[1]; k++)
				{
					if (population[index[i]].chorm[j] == population[index[i]].chorm[k])
					{
						population[index[i]].chorm[j] = population[index[i + 1]].chorm[k];
						k = point[0]-1;
					}
				}
			}
			for (j = point[1]; j < LENGTH; j++)
			{
				for (k = point[0]; k < point[1]; k++)
				{
					if (population[index[i]].chorm[j] == population[index[i]].chorm[k])
					{
						population[index[i]].chorm[j] = population[index[i + 1]].chorm[k];
						k = point[0]-1;
					}
				}
			}
			for (j = 0; j < point[0]; j++)
			{
				for (k = point[0]; k < point[1]; k++)
				{
					if (population[index[i+1]].chorm[j] == population[index[i + 1]].chorm[k])
					{
						population[index[i + 1]].chorm[j] = population[index[i]].chorm[k];
						k = point[0]-1;
					}
				}
			}
			for (j = point[1]; j < LENGTH; j++)
			{
				for (k = point[0]; k < point[1]; k++)
				{
					if (population[index[i + 1]].chorm[j] == population[index[i + 1]].chorm[k])
					{
						population[index[i + 1]].chorm[j] = population[index[i]].chorm[k];
						k = point[0]-1;
					}
				}
			}
		}
	}
}


void MutationOperation(void)
{

	int i, j, k, m, n;
	int DDT[LENGTH][LENGTH], LAT[LENGTH][LENGTH], BCT[LENGTH][LENGTH];
	int detain, detaout, Sin, alpha, beta;
	int sum1, sum2;
	int inverseSbox[LENGTH];
	double p;
	int position[6], back[6],mark[LENGTH];
	struct individual minf;


	for (i = 0; i < POPSIZE; i++)
	{
		
		p = rand() / (double)(RAND_MAX);

		if (p < Pm)
		{

			for (j = 0; j < LENGTH; j++)
			{
				for (k = 0; k < LENGTH; k++)
				{
					DDT[j][k] = 0;
					BCT[j][k] = 0;
					LAT[j][k] = -pow((double)2, BITS - 1);
				}
			}

			for (detain = 0; detain < LENGTH; detain++)
			{
				for (Sin = 0; Sin < LENGTH; Sin++)
				{
					detaout = population[i].chorm[Sin] ^ population[i].chorm[Sin^detain];
					DDT[detain][detaout] ++;
				}
			}

			for (alpha = 0; alpha < LENGTH; alpha++)
			{
				for (beta = 0; beta < LENGTH; beta++)
				{
					for (Sin = 0; Sin < LENGTH; Sin++)
					{
						if (innerproduct(Sin, alpha) == innerproduct(population[i].chorm[Sin], beta))
							LAT[alpha][beta] ++;
					}
				}
			}

			for (j = 0; j<LENGTH; j++)
			{
				inverseSbox[population[i].chorm[j]] = j;
			}
			for (detain = 0; detain < LENGTH; detain++)
			{
				for (detaout = 0; detaout< LENGTH; detaout++)
				{
					for (Sin = 0; Sin < LENGTH; Sin++)
					{
						if (detain == ((inverseSbox[population[i].chorm[Sin] ^ detaout]) ^ (inverseSbox[population[i].chorm[Sin^detain] ^ detaout])))
						{
							BCT[detain][detaout] ++;
						}
					}
				}
			}

			population[i].duniformity = 0;
			for (j = 0; j < LENGTH; j++)
			{
				for (k = 0; k < LENGTH; k++)
				{
					if ((DDT[j][k] > population[i].duniformity) && ((j + k) != 0))
						population[i].duniformity = DDT[j][k];
				}
			}
			population[i].linearity = 0;
			for (j = 0; j < LENGTH; j++)
			{
				for (k = 0; k < LENGTH; k++)
				{
					if ((abs(LAT[j][k]) > abs(population[i].linearity)) && ((j + k) != 0))
						population[i].linearity = LAT[j][k];
				}
			}
			population[i].linearity = 2 * population[i].linearity;
			population[i].buniformity = 0;
			for (j = 0; j < LENGTH; j++)
			{
				for (k = 0; k < LENGTH; k++)
				{
					if ((BCT[j][k] > population[i].buniformity) && (j != 0) && (k != 0))
						population[i].buniformity = BCT[j][k];
				}
			}

			population[i].fitness = population[i].duniformity + abs(population[i].linearity) + population[i].buniformity;




			minf = population[i];
			for (j = 0; j < LENGTH; j++)
			{
				mark[j] = 0;
			}
			for (j = 0; j < 6; j++)
			{
				position[j] = rand() % LENGTH;
				if (mark[position[j]] == 0)
				{
					mark[position[j]] = 1;
				}
				else
				{
					j--;
				}
			}
			
			for (j = 0; j < 6; j++)
			{
				back[j] = population[i].chorm[position[j] ];
			}
			for (m = 1; m < 720; m++)
			{
				for (j = 0; j < 6; j++)
				{
					population[i].chorm[position[j]] = back[C[m][j]];
				}

				for (j = 0; j < LENGTH; j++)
				{
					for (k = 0; k < LENGTH; k++)
					{
						DDT[j][k] = 0;
						BCT[j][k] = 0;
						LAT[j][k] = -pow((double)2, BITS - 1);
					}
				}

				for (detain = 0; detain < LENGTH; detain++)
				{
					for (Sin = 0; Sin < LENGTH; Sin++)
					{
						detaout = population[i].chorm[Sin] ^ population[i].chorm[Sin^detain];
						DDT[detain][detaout] ++;
					}
				}

				for (alpha = 0; alpha < LENGTH; alpha++)
				{
					for (beta = 0; beta < LENGTH; beta++)
					{
						for (Sin = 0; Sin < LENGTH; Sin++)
						{
							if (innerproduct(Sin, alpha) == innerproduct(population[i].chorm[Sin], beta))
								LAT[alpha][beta] ++;
						}
					}
				}

				for (j = 0; j<LENGTH; j++)
				{
					inverseSbox[population[i].chorm[j]] = j;
				}
				for (detain = 0; detain < LENGTH; detain++)
				{
					for (detaout = 0; detaout< LENGTH; detaout++)
					{
						for (Sin = 0; Sin < LENGTH; Sin++)
						{
							if (detain == ((inverseSbox[population[i].chorm[Sin] ^ detaout]) ^ (inverseSbox[population[i].chorm[Sin^detain] ^ detaout])))
							{
								BCT[detain][detaout] ++;
							}
						}
					}
				}

				population[i].duniformity = 0;
				for (j = 0; j < LENGTH; j++)
				{
					for (k = 0; k < LENGTH; k++)
					{
						if ((DDT[j][k] > population[i].duniformity) && ((j + k) != 0))
							population[i].duniformity = DDT[j][k];
					}
				}
				population[i].linearity = 0;
				for (j = 0; j < LENGTH; j++)
				{
					for (k = 0; k < LENGTH; k++)
					{
						if ((abs(LAT[j][k]) > abs(population[i].linearity)) && ((j + k) != 0))
							population[i].linearity = LAT[j][k];
					}
				}
				population[i].linearity = 2 * population[i].linearity;
				population[i].buniformity = 0;
				for (j = 0; j < LENGTH; j++)
				{
					for (k = 0; k < LENGTH; k++)
					{
						if ((BCT[j][k] > population[i].buniformity) && (j != 0) && (k != 0))
							population[i].buniformity = BCT[j][k];
					}
				}

				population[i].fitness = population[i].duniformity +  abs(population[i].linearity) + population[i].buniformity;

				if (population[i].fitness < minf.fitness)
				{
					minf = population[i];
				}

			}

			population[i] = minf;

		}

	}

}