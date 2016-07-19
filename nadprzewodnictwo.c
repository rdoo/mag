#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// STALE
#define nm2au 18.89726133921252
#define eV2au 0.03674932587122423
#define k2au 0.00316681520371153
#define mili 0.001

const double masa_e = 1.0; // masa elektronu w j. a.
const double h = 1.0; // h kreslone w j. a.
const double stala_kb = 3.16681520371153e-6; // stala boltzmanna w j. a.

// PARAMETRY PROGRAMU
const int N = 20; // liczba pasm
const double prog_samouzgodnienia = 1.0e-8 * eV2au; // jezeli roznica pomiedzy poprzednia delta i nastepna jest mniejsza to program stopuje
const double delta0 = 0.25 * mili * eV2au; // poczatkowa przerwa nadprzewodzaca
double T = 0.1; // temperatura
const double potencjal_chem_09eV = 0.9 * eV2au; // poczatkowy potencjal chemiczny
const double gN0 = 0.18; // potrzebne do stalej oddzialywania g
const double EDebye = 32.31 * mili * eV2au; // energia Debye'a
const double L_min = nm2au; // poczatkowa grubosc warstwy
const double L_max = 5 * nm2au; // koncowa grubosc warstwy
const double dL = nm2au / 20;
double potencjal_chem = 0.9 * eV2au; // potencjal chemiczny
double L; // aktualna grubosc warstwy (zainicjalizowana w glownej petli)
double g; // stala oddzialywania elektron-fonon (zainicjalizowana w mainie)
double kF; // wektor Fermiego

const int niejednorodnosc = 1; // czy brac pod uwage niejednorodnosc powierzchni

// PARAMETRY CALKOWANIA
const double dz = 0.01 * nm2au; // przedzial calkowania 

const double k_min = 0.0 / nm2au; // k min do calkowania po k
const double k_max = 8.0 / nm2au; // k max do calkowania po k
const double dk = 0.001 / nm2au; // dk do calkowania po k


const double T_min = 0.1;
const double T_max = 10;
const double dT = 0.1;

const int max_liczba_iteracji = 1e6;
const double prog_akceptacji_Tc = 1e-4;


// pomocnicza funkcja wypisujaca na ekran
void wypisz(char* string, double value) {
	if (fabs(value) > 0.001 && fabs(value) < 1000) {
		printf("%s: %.3f\n", string, value);
	} else {
		printf("%s: %.2e\n", string, value);
	}
}

// pomocnicza funkcja sumujaca tablice
double sumaTablicy(double* tablica) {
	double suma = 0.;
	int i;
	for (i = 0; i < N; i++) {
		suma += tablica[i];
	}
	return suma;
}

double wartoscFunkcjiFalowejStudni(double x, int i) {
	return sqrt(2. / L) * sin(i * M_PI * x / L);
}

// calkowanie metoda prostokatow
double wartoscPoczatkowejDelty(int i) {
	double zp = 0.; // z poczatkowe
	double zk = L; // z koncowe
	
	double calka = 0.;
	double z;
	for (z = zp; z <= zk; z += dz) {
		calka += wartoscFunkcjiFalowejStudni(z, i) * delta0 * wartoscFunkcjiFalowejStudni(z, i);
	}
	calka *= dz;
	return calka;
}

double wartoscC(int i, int j) {
	double zp = 0.; // z poczatkowe
	double zk = L; // z koncowe
	
	double calka = 0.;
	double z;
	for (z = zp; z <= zk; z += dz) {
		calka += wartoscFunkcjiFalowejStudni(z, j) * wartoscFunkcjiFalowejStudni(z, j) * wartoscFunkcjiFalowejStudni(z, i) * wartoscFunkcjiFalowejStudni(z, i);
	}
	calka *= dz;
	return calka;
}

double wartoscKsi(int i) {
	return i * i * M_PI * M_PI / (2 * masa_e * L * L);
}

double wartoscEnergii(double k, double deltai, double ksii) {
	double czynnik1 = ksii + k * k / 2 - potencjal_chem; // 1. czynnik sumy pod pierwiastkiem
	return sqrt(czynnik1 * czynnik1 + deltai * deltai);
}

double wartoscRozkladuFermiego(double E) {
	return 1.0 / (1.0 + exp(E / (stala_kb * T)));
}

double wartoscDelta(double* poprzednia_delta, double (*C)[N], int j, double Ekin_add) {
	double kp = k_min; // k poczatkowe
	double kk = k_max; // k koncowe

	double calka = 0.;

	int i;
	double k, E, Ekin;
	for (k = kp; k <= kk; k += dk) {
		
		for (i = 0; i < N; i++) {
			Ekin = (i+1)*(i+1)*M_PI*M_PI/(2*masa_e*L*L)+k*k/2.0/masa_e-potencjal_chem;

			if (niejednorodnosc) {
				Ekin += Ekin_add;
			}

			E=sqrt(Ekin*Ekin+poprzednia_delta[i]*poprzednia_delta[i]);
			//printf("%e\n",Ekin/eV2au);
			if (fabs(Ekin) <= EDebye) {				
				calka += (g/(2.0*M_PI))*C[i][j]*poprzednia_delta[i]/(2.0*E)*(1.0 - 2.0 * wartoscRozkladuFermiego(E))*k*dk;
			}
		}
		
	}
	return calka;
}


// sprawdza roznice pomiedzy sumami delt poprzedniej i nastepnej
int sprawdzRoznice(double* poprzednia_delta, double* nastepna_delta) {
	double poprzednia_suma = sumaTablicy(poprzednia_delta);	
	double nastepna_suma = sumaTablicy(nastepna_delta);

	//wypisz("poprzednia suma wynosi", poprzedniaSuma);
	//wypisz("nastepna suma wynosi", nastepnaSuma);
	
	// roznica pomiedzy suma poprzedniej i nastepnej delty
	// TODO jakos inaczej to sprawdzac??
	double roznica = fabs(nastepna_suma - poprzednia_suma);
	
	if (roznica > prog_samouzgodnienia) {
		//wypisz("Roznica jest duza i wynosi", roznica);
		return 1;
	}
	//wypisz("Roznica jest mala i wynosi", roznica);
	return 0;
}

double gestoscElektronow(double potencjal_chem, double* poczatkowa_delta) {
	double calka = 0.0;

	double k, z;
	int i;
	for (k = k_min; k <= k_max; k += dk) {
		double sumaN = 0.0;
		for (i = 0; i < N; i++) {
			double calkaZ = 0.0;
			for (z = 0; z <= L; z += dz) {
				double energia_kinetyczna = wartoscKsi(i+1) + h*h*k*k/(2*masa_e) - potencjal_chem;
				double energia_calkowita = sqrt(energia_kinetyczna * energia_kinetyczna + poczatkowa_delta[i] * poczatkowa_delta[i]);
				double suma = energia_kinetyczna + energia_calkowita;
				double pierwiastek = sqrt(suma * suma + poczatkowa_delta[i] * poczatkowa_delta[i]);
				double wartosc_f_falowej = wartoscFunkcjiFalowejStudni(z, i+1);

				double U = suma / pierwiastek;
				double V = poczatkowa_delta[i] / pierwiastek;
				double fE = wartoscRozkladuFermiego(energia_calkowita);

				calkaZ += U*U*wartosc_f_falowej*wartosc_f_falowej*fE+V*V*wartosc_f_falowej*wartosc_f_falowej*(1-fE);
			}
			calkaZ *= dz;

			sumaN += calkaZ;
		}
		calka += sumaN * k;
	}
	return calka * dk / M_PI / L;
}

double potencjalChemicznyMetodaBisekcji(double poczatkowa_gestosc_elektronow, double* poczatkowa_delta) {
	double a = 0.0 * eV2au;
	double b = 2.0 * eV2au;
	double prog_dokladnosci = 1e-3 * eV2au;
	
	while (fabs(a - b) > prog_dokladnosci) {
		double x = (a + b) / 2;
		double gestosc1 = gestoscElektronow(x, poczatkowa_delta);
		double gestosc2 = gestoscElektronow(a, poczatkowa_delta);
		if ((gestosc1 - poczatkowa_gestosc_elektronow) * (gestosc2 - poczatkowa_gestosc_elektronow) < 0) {
			b = x;
		} else {
			a = x;
		}
	}

	return (a + b) / 2;
}

double gestoscElektronowWModelu3D() {
	double prog_dokladnosci = 1e-6;

	double stala_przed_calka = 1.0 / 2.0 / M_PI / M_PI * pow(2.0 * masa_e / h / h * stala_kb * T, 3.0 / 2.0);
	double stala_w_calce = potencjal_chem_09eV / stala_kb / T;

	double calka = 0.0;

	double x = 0.0;
	double dx = dz;
	double przyczynek;
	do {
		x += dx;
		przyczynek = sqrt(x) / (1 + exp(x - stala_w_calce));
		calka += przyczynek;
	} while (przyczynek > prog_dokladnosci);

	return stala_przed_calka * calka * dx;
}


void obliczanieDyspersji(double* delta_nadprzewodzaca) {
	FILE *plik_energia_od_k;
	char nazwa_pliku[64];
	snprintf(nazwa_pliku, sizeof nazwa_pliku, "dane/dyspersja_dla_L_%.2f.txt", L/nm2au);
	plik_energia_od_k = fopen(nazwa_pliku, "w");

	double energia_kinetyczna[N];
	double energia_calkowita[N];

	double k;
	for (k = k_min; k <= k_max; k += dk) {
		fprintf(plik_energia_od_k, "%.3f", k*nm2au);

		int i;
		for (i = 0; i < 4; i++) {
			energia_kinetyczna[i] = wartoscKsi(i+1) + h*h*k * k / (2 * masa_e) - potencjal_chem;
			energia_calkowita[i] = sqrt(energia_kinetyczna[i] * energia_kinetyczna[i] + delta_nadprzewodzaca[i] * delta_nadprzewodzaca[i]);
		}

		for (i = 0; i < 4; i++) {
			fprintf(plik_energia_od_k, " %.20f", energia_calkowita[i]/eV2au/mili);
		}

		for (i = 0; i < 4; i++) {
			fprintf(plik_energia_od_k, " %.20f", energia_kinetyczna[i]/eV2au/mili);
		}

		fprintf(plik_energia_od_k, "\n");
	}

	fclose(plik_energia_od_k);
}

void obliczanieRozkladuDeltaOdZ(double* delta_nadprzewodzaca) {
	FILE *plik_delta_od_z;
	char nazwa_pliku[64];
	snprintf(nazwa_pliku, sizeof nazwa_pliku, "dane/delta_od_z_dla_L_%.2f.txt", L/nm2au);
	plik_delta_od_z = fopen(nazwa_pliku, "w");

	double kp = k_min; // k poczatkowe
	double kk = k_max; // k koncowe


	int i;
	double z, k, E, Ekin;
	for (z = 0.0; z <= L; z += dz) {
		double calka = 0.;

		for (k = kp; k <= kk; k += dk) {
			
			for (i = 0; i < N; i++) {
				Ekin = (i+1)*(i+1)*M_PI*M_PI/(2*masa_e*L*L)+k*k/2.0/masa_e-potencjal_chem;
				E=sqrt(Ekin*Ekin+delta_nadprzewodzaca[i]*delta_nadprzewodzaca[i]);
				//printf("%e\n",Ekin/eV2au);
				if (fabs(Ekin) <= EDebye) {				
					calka += (g/(2.0*M_PI))*delta_nadprzewodzaca[i]/(2.0*E)*wartoscFunkcjiFalowejStudni(z, i+1)*wartoscFunkcjiFalowejStudni(z, i+1)*(1.0 - wartoscRozkladuFermiego(E))*k*dk;
				}
			}
		}
		fprintf(plik_delta_od_z, "%.3f %.20f\n", z/nm2au, calka/eV2au/mili);
	}

	fclose(plik_delta_od_z);
}

int main() {
	srand(time(NULL));

	kF=sqrt(2.0*masa_e*potencjal_chem);
	g = gN0 / (masa_e * kF / (2. * M_PI * M_PI));

	T = T_min;

	// wypisz("Prog samouzgodnienia wynosi", progSamouzgodnienia);
	// wypisz("Delta0 wynosi", delta0);
	// wypisz("EDebye wynosi", EDebye);
	// wypisz("g elektron-fonon", g);
	// wypisz("wektor Fermiego", kF);
	
	FILE *plik_delta_od_L;
	plik_delta_od_L = fopen("dane/delta_od_L.txt", "w");

	FILE *plik_potencjal_od_L;
	plik_potencjal_od_L = fopen("dane/potencjal_od_L.txt", "w");

	FILE *plik_Tc_od_L;
	plik_Tc_od_L = fopen("dane/Tc_od_L.txt", "a");

	int licznik_petli; // jesli chcemy cale obiczenia wykonac kilkukrotnie np. do obliczenia sredniej
	for (licznik_petli = 0; licznik_petli < 10; licznik_petli++) {
		wypisz("Numer petli", licznik_petli);

		double poczatkowa_gestosc_elektronow = gestoscElektronowWModelu3D();
		wypisz("Poczatkowa gestosc elektronow na cm^3", poczatkowa_gestosc_elektronow*nm2au*nm2au*nm2au*1e21);
		
		// glowna petla liczaca delty dla roznych grubosci L
		for (L = L_min; L <= L_max; L += dL) {
		  

		  	T = T_min;
		  
			//L=2*nm2au;
			wypisz("Obliczenia dla L w nm", L/nm2au);
			int i, j;


			// tablica przechowujaca wartosc poczatkowa DELTAi
			double poczatkowa_delta[N];
			for (i = 0; i < N; i++) {
				poczatkowa_delta[i] = wartoscPoczatkowejDelty(i + 1);
				//wypisz("Poczatkowa delta", wartoscPoczatkowejDelty(i + 1));
			}
			
			
			potencjal_chem = potencjalChemicznyMetodaBisekcji(poczatkowa_gestosc_elektronow, poczatkowa_delta); // potencjal czyli mi

			wypisz("Potencjal chemiczny w eV", potencjal_chem/eV2au);

			
			fprintf(plik_potencjal_od_L, "%.2f %.20f\n", L/nm2au, potencjal_chem/eV2au);
			
			


			// tablica przechowujaca wartosc stalych Cij
			double C[N][N];
			for (i = 0; i < N; i++) {
				for (j = 0; j < N; j++) {
						C[i][j] = wartoscC(i + 1, j + 1);
						//wypisz("C", wartoscC(i + 1, j + 1));
				}
			}

			double Ekin_add = 0.;

			if (niejednorodnosc) {
				double alfa = (double)rand() / (double)RAND_MAX * 2 - 1; // liczba losowa z przedzialu -1 do 1 TODO: dobrze?
				Ekin_add = alfa * (i+1)*(i+1)*M_PI*M_PI/(masa_e*L*L*L) * 0.286 * nm2au; // 0.286 to 1 ML dla Pb
			}


			FILE *plik_delta_od_T;
			char nazwa_pliku[64];
			snprintf(nazwa_pliku, sizeof nazwa_pliku, "dane/delta_od_T_dla_L_%.2f.txt", L/nm2au);
			plik_delta_od_T = fopen(nazwa_pliku, "w");

			int pierwsza_petla = 1;
			for (T = T_min; T <= T_max; T += dT) {
				wypisz("Aktualna temperatura", T);

				// ALGORYTM SAMOUZGODNIENIA TODO nie wiem czy poprawny
				double poprzednia_delta[N];
				double poprzednia_delta_kopia[N];
				double nastepna_delta[N];
				
				for (i = 0; i < N; i++) {
					poprzednia_delta[i] = poczatkowa_delta[i];
				}

				int liczba_iteracji = 0;
				
				do {
					liczba_iteracji++;
					//wypisz("liczba_iteracji", liczba_iteracji);

					if (liczba_iteracji > max_liczba_iteracji) {
						break;
					}

					for (i = 0; i < N; i++) {
						poprzednia_delta_kopia[i] = poprzednia_delta[i];
					}
					
					//for (i = 0; i < N; i++) {
					//	printf("%e ", poprzedniaDelta[i]);
					//}
					//printf("\n");
					
					for (i = 0; i < N; i++) {
						nastepna_delta[i] = wartoscDelta(poprzednia_delta, C, i, Ekin_add);
					}
					
					for (i = 0; i < N; i++) {
						poprzednia_delta[i] = nastepna_delta[i];
					}
					
				} while (sprawdzRoznice(poprzednia_delta_kopia, nastepna_delta));

				//double sumaKoncowaDelty = sumaTablicy(nastepnaDelta);
				//wypisz("Koncowy wynik delty", sumaKoncowaDelty);
				wypisz("Liczba iteracji samouzgodnienia", liczba_iteracji);
				fprintf(plik_delta_od_T, "%.2f %.20f\n", T, nastepna_delta[0]/eV2au/mili);

				if (pierwsza_petla == 1 && !niejednorodnosc) { // nie licz tych rzeczy jesli liczymy niejednorodnosc

					fprintf(plik_delta_od_L, "%.2f %.20f\n", L/nm2au, nastepna_delta[0]/eV2au/mili);

					obliczanieDyspersji(nastepna_delta);
					obliczanieRozkladuDeltaOdZ(nastepna_delta);
					pierwsza_petla = 0;
				}

				if (nastepna_delta[0]/eV2au/mili < prog_akceptacji_Tc || liczba_iteracji > max_liczba_iteracji) {
					fprintf(plik_Tc_od_L, "%.2f %.20f\n", L/nm2au, T);
					break;
				}

				if (niejednorodnosc) { // jesli liczymy niejednorodnosc to nie interesuja nas pozostale temperatury
					break;
				}

			}
			fclose(plik_delta_od_T);


			if (!niejednorodnosc) { // przerwij petle jesli nie jest liczona niejednorosc
				break;
			}
		}

	}

	fclose(plik_delta_od_L);
	fclose(plik_potencjal_od_L);
	fclose(plik_Tc_od_L);

	system("gnuplot plot_delta");
	system("gnuplot plot_dyspersja");
	system("gnuplot plot_potencjal");
	system("gnuplot plot_delta_od_z");
	system("gnuplot plot_delta_od_T");
	system("gnuplot plot_Tc_od_L");
}

