#include <iostream>
#include <fstream>
using namespace std;

const int mx = 201, my = 121, ma = max(mx, my), mi = 1;
int itk, ish, isht, it; //итерации и шаг

double tet[mx+1][my + 1], tets[mx + 1][my + 1], tet1[mx + 1][my + 1], tet2[mx + 1][my + 1];
double vi[mx + 1][my + 1], vis[mx + 1][my + 1], vi1[mx + 1][my + 1], vi2[mx + 1][my + 1];
double ft[mx + 1][my + 1], fts[mx + 1][my + 1], ft1[mx + 1][my + 1], ft2[mx + 1][my + 1];
double Ux[mx + 1][my + 1], Uy[mx + 1][my + 1];
double a[ma + 1], b[ma + 1], c[ma + 1], d[ma + 1], e[ma + 1];

double Re, Pr, dt, dl_trub, dx, dy, sumft, sumvi, sumtet, Gr, x, y;
string sim, name;

void Tom(int i0, int iN, double(&a)[ma + 1], double(&b)[ma + 1], double(&c)[ma + 1], double(&d)[ma + 1], double(&e)[ma + 1]) {
	double alf[ma + 1], bet[ma + 1];

	alf[i0] = c[i0] / b[i0];
	bet[i0] = d[i0] / b[i0];

	for (int i = i0 + 1; i <= iN - 1; i++) {
		alf[i] = c[i] / (b[i] - a[i] * alf[i - 1]);
	}

	for (int i = i0 + 1; i <= iN; i++) {
		bet[i] = (d[i] + a[i] * bet[i - 1]) / (b[i] - a[i] * alf[i - 1]);
	}

	e[iN] = bet[iN];

	for (int i = iN - 1; i >= i0; i--) {
		e[i] = alf[i] * e[i + 1] + bet[i];
	}
}


int main()
{
	system("chcp 1251 > 0");

	dt = 0.001;
	Gr = 0.0;
	Re = 10.0;
	Pr = 1.0;
	double x1 = 0.5, y1 = 0.4, y2 = 0.6;
	double UxVerh = 1.0;
	double dl_trub = 4.0;

	dx = dl_trub / (mx - 1);
	dy = 1.0 / (my - 1);
	double dx2 = dx * dx;
	double dy2 = dy * dy;

	int i1 = static_cast<int>(round(x1 / dx) + 1);
	int j1 = static_cast<int>(round(y1 / dy) + 1);
	int j2 = static_cast<int>(round(y2 / dy) + 1);

	x1 = (i1 - 1) * dx;
	y1 = (j1 - 1) * dy;
	y2 = (j2 - 1) * dy;

	//ГРАНИЧНЫЕ УСЛОВИЯ
		//Температура
	for (int j = 2; j <= j1 - 1; j++) {
		tet[1][j] = 0.0;
	}
	for (int j = j2+1; j <= my - 1; j++) {
		tet[1][j] = 1.0;
	}

	for (int i = 1; i <= mx; i++) {
		for (int j = 1; j <= my; j++) {
			tets[i][j] = tet[i][j];
		}
	}

	//Скорость
	for (int j = 2; j <= j1 - 1; j++) {
		Ux[1][j] = 1.0;
	}
	for (int j = j2+1; j <= my - 1; j++) {
		Ux[1][j] = UxVerh;
	}

	//Функция тока
	for (int j = 1; j <= j1; j++) {
		y = (j-1)*dy;
		ft[1][j] = y;
	}

	for (int i = 2; i <= i1; i++) {
		ft[i][j1] = ft[1][j1];
	}

	for (int j = j1; j <= j2; j++) {
		ft[i1][j] = ft[1][j1];
	}

	for (int i = 1; i <= i1; i++) {
		ft[i][j2] = ft[1][j1];
	}

	for (int j = j2+1; j <= my; j++) {
		y = (j-1)*dy;
		ft[1][j] = ft[1][j1] + (y - y2) * UxVerh;
	}

	for (int i = 2; i <= mx; i++) {
		ft[i][my] = ft[1][my];
		ft[i][1] = 0.0;
	}

	for (int i = 1; i <= mx; i++) {
		for (int j = 1; j <= my; j++) {
			fts[i][j] = ft[i][j];
		}
	}
	//Расчет
	cout << "Введите количество итераций: ";
	cin >> itk;

	cout << "Введите шаг: ";
	cin >> ish;
	isht = ish;
	it = 0;

	while (it <= itk) {
		it++;

		//Функция тока
		{
			//Прогонка по Х
			{
				//Область Ox_1
				for (int j = 2; j <= j1 - 1; j++) {
					for (int i = 2; i <= mx - 1; i++) {
						a[i] = 1. / dx2;
						c[i] = 1. / dx2;
						b[i] = 1. / dt + 2. / dx2;
						d[i] = fts[i][j] / dt - vi[i][j];
					}

					a[1] = 0.0;       b[1] = 1.0;       c[1] = 0.0;       d[1] = fts[1][j];
					a[mx] = 1.0;	   b[mx] = 1.0;		c[mx] = 0.0;	  d[mx] = 0.0;

					Tom(1, mx, a, b, c, d, e);

					for (int i = 1; i <= mx; i++)
						ft1[i][j] = e[i];
				}

				//Область Ox_2
				for (int j = j1; j <= j2; j++) {
					for (int i = i1 + 1; i <= mx - 1; i++) {
						a[i] = 1. / dx2;
						c[i] = 1. / dx2;
						b[i] = 1. / dt + 2. / dx2;
						d[i] = fts[i][j] / dt - vi[i][j];
					}

					a[i1] = 0.0;  b[i1] = 1.0;       c[i1] = 0.0;       d[i1] = fts[i1][j];
					a[mx] = 1.0;  b[mx] = 1.0;       c[mx] = 0.0;       d[mx] = 0.0;

					Tom(i1, mx, a, b, c, d, e);

					for (int i = i1; i <= mx; i++)
						ft1[i][j] = e[i];
				}

				//Область Ox_3
				for (int j = j2+1; j <= my - 1; j++) {
					for (int i = 2; i <= mx - 1; i++) {
						a[i] = 1. / dx2;
						c[i] = 1. / dx2;
						b[i] = 1. / dt + 2. / dx2;
						d[i] = fts[i][j] / dt - vi[i][j];
					}

					a[1] = 0.0;       b[1] = 1.0;       c[1] = 0.0;       d[1] = fts[1][j];
					a[mx] = 1.0;	  b[mx] = 1.0;		c[mx] = 0.0;	  d[mx] = 0.0;

					Tom(1, mx, a, b, c, d, e);

					for (int i = 1; i <= mx; i++)
						ft1[i][j] = e[i];
				}
			}

			//Прогонка по Y
			{
				//Область Oy_1
				for (int i = 2; i <= i1; i++) {
					for (int j = 2; j <= j1 - 1; j++) {
						a[j] = 1. / dy2;
						c[j] = 1. / dy2;
						b[j] = 1. / dt + 2. / dy2;
						d[j] = ft1[i][j] / dt;
					}

					a[1] = 0.0;       b[1] = 1.0;       c[1] = 0.0;       d[1] = 0.0;
					a[j1] = 0.0;  b[j1] = 1.0;  c[j1] = 0.0;  d[j1] = fts[i][j1];

					Tom(1, j1, a, b, c, d, e);

					for (int j = 1; j <= j1; j++)
						ft[i][j] = e[j];
				}

				//Область Oy_2
				for (int i = 2; i <= i1; i++)
				{
					for (int j = j2+1; j <= my - 1; j++) {
						a[j] = 1. / dy2;
						c[j] = 1. / dy2;
						b[j] = 1. / dt + 2. / dy2;
						d[j] = ft1[i][j] / dt;
					}


					a[j2] = 0.0;     b[j2] = 1.0;    c[j2] = 0.0;   d[j2] = fts[i][j2];
					a[my] = 0.0;     b[my] = 1.0;    c[my] = 0.0;   d[my] = fts[i][my];

					Tom(j2, my, a, b, c, d, e);

					for (int j = j2; j <= my; j++)
						ft[i][j] = e[j];
				}

				//Область Oy_3
				for (int i = i1+1; i <= mx - 1; i++)
				{
					for (int j = 2; j <= my - 1; j++) {
						a[j] = 1. / dy2;
						c[j] = 1. / dy2;
						b[j] = 1. / dt + 2. / dy2;
						d[j] = ft1[i][j] / dt;
					}


					a[1] = 0.0;       b[1] = 1.0;       c[1] = 0.0;       d[1] = fts[i][1];
					a[my] = 0.0;	  b[my] = 1.0;		c[my] = 0.0;	  d[my] = fts[i][my];

					Tom(1, my, a, b, c, d, e);

					for (int j = 1; j <= my; j++)
						ft[i][j] = e[j];
				}
			}
			//Границы
			for (int j = 2; j <= my - 1; j++) {
				ft[mx][j] = ft[mx - 1][j];
			}

			for (int i = 1; i <= mx; i++) {
				for (int j = 1; j <= my; j++) {
					ft2[i][j] = ft[i][j] - fts[i][j];
				}
			}
		}

		//СКОРОСТЬ
		{
			for (int j = 2; j <= my - 1; j++) {
				for (int i = 2; i <= mx - 1; i++) {
					if (!(i <= i1 && j >= j1 && j <= j2)) {
						Ux[i][j] = (ft[i][j + 1] - ft[i][j - 1]) / 2. / dy;
						Uy[i][j] = -(ft[i + 1][j] -	ft[i - 1][j]) / 2. / dx;
					}
				}
			}

			for (int j = 2; j <= my - 1; j++) {
				Ux[mx][j] = Ux[mx - 1][j];
				Uy[mx][j] = Uy[mx - 1][j];
			}

			for (int i = 1; i <= i1; i++) {
				Ux[i][j1] = 0.;
				Uy[i][j1] = 0.;
				Ux[i][j2] = 0.;
				Uy[i][j2] = 0.;
			}

			for (int j = j1; j <= j2; j++) {
				Ux[i1][j] = 0.;
				Uy[i1][j] = 0.;
			}

			for (int j = j2+1; j <= my - 1; j++) {
				Uy[1][j] = Uy[2][j];
			}

			for (int j = 2; j <= j1 - 1; j++) {
				Uy[1][j] = Uy[2][j];
			}
		}

		//Вихрь
		{
			//Прогонка по Х
			{
				//Область Ox_1
				for (int j = 2; j <= j1 - 1; j++) {
					for (int i = 2; i <= mx - 1; i++) {
						double aUx = abs(Ux[i][j]);
						a[i] = (aUx + Ux[i][j]) / 2. / dx + 1. / Re / dx2;
						c[i] = (aUx - Ux[i][j]) / 2. / dx + 1. / Re / dx2;
						b[i] = 1. / dt + 2. / Re / dx2 + aUx / dx;
						d[i] = vis[i][j] / dt - Gr / Re / Re * (tet[i][j + 1] - tet[i][j - 1]) / 2. / dy;
					}

					a[1] = 0.0;     b[1] = 1.0;       c[1] = 0.0;       d[1] = 0.0;
					a[mx] = 1.0;	b[mx] = 1.0;	  c[mx] = 0.0;      d[mx] = 0.0;

					Tom(1, mx, a, b, c, d, e);

					for (int i = 1; i <= mx; i++)
						vi1[i][j] = e[i];
				}

				//Область Ox_2
				for (int j = j1; j <= j2; j++) {
					for (int i = i1 + 1; i <= mx - 1; i++) {
						double aUx = abs(Ux[i][j]);
						a[i] = (aUx + Ux[i][j]) / 2. / dx + 1. / Re / dx2;
						c[i] = (aUx - Ux[i][j]) / 2. / dx + 1. / Re / dx2;
						b[i] = 1. / dt + 2. / Re / dx2 + aUx / dx;
						d[i] = vis[i][j] / dt - Gr / Re / Re * (tet[i][j + 1] - tet[i][j - 1]) / 2. / dy;
					}

					a[i1] = 0.0;  b[i1] = 1.0;       c[i1] = 0.0;       d[i1] = 2 * (ft[i1 + 1][j] - ft[i1][j]) / dx2;
					a[mx] = 1.0;  b[mx] = 1.0;       c[mx] = 0.0;       d[mx] = 0.0;

					Tom(i1, mx, a, b, c, d, e);

					for (int i = i1; i <= mx; i++)
						vi1[i][j] = e[i];
				}

				//Область Ox_3
				for (int j = j2 + 1; j <= my - 1; j++) {
					for (int i = 2; i <= mx - 1; i++) {
						double aUx = abs(Ux[i][j]);
						a[i] = (aUx + Ux[i][j]) / 2. / dx + 1. / Re / dx2;
						c[i] = (aUx - Ux[i][j]) / 2. / dx + 1. / Re / dx2;
						b[i] = 1. / dt + 2. / Re / dx2 + aUx / dx;
						d[i] = vis[i][j] / dt - Gr / Re / Re * (tet[i][j + 1] - tet[i][j - 1]) / 2. / dy;
					}

					a[1] = 0.0;       b[1] = 1.0;       c[1] = 0.0;       d[1] = 0.0;
					a[mx] = 1.0;	  b[mx] = 1.0;      c[mx] = 0.0;      d[mx] = 0.0;

					Tom(1, mx, a, b, c, d, e);

					for (int i = 1; i <= mx; i++)
						vi1[i][j] = e[i];
				}
			}

			//Прогонка по Y
			{
				//Область Oy_1
				for (int i = 2; i <= i1; i++)
				{
					for (int j = 2; j <= j1 - 1; j++) {
						double aUy = abs(Uy[i][j]);
						a[j] = (aUy + Uy[i][j]) / 2. / dy + 1. / Re / dy2;
						c[j] = (aUy - Uy[i][j]) / 2. / dy + 1. / Re / dy2;
						b[j] = 1. / dt + 2. / Re / dy2 + aUy / dy;
						d[j] = vi1[i][j] / dt;
					}

					a[1] = 0.0;      b[1] = 1.0;     c[1] = 0.0;     d[1] = 2 * (ft[i][2] - ft[i][1]) / dy2;
					a[j1] = 0.0;	 b[j1] = 1.0;	 c[j1] = 0.0;	 d[j1] = 2 * (ft[i][j1 - 1] - ft[i][j1]) / dy2;

					Tom(1, j1, a, b, c, d, e);

					for (int j = 1; j <= j1; j++)
						vi[i][j] = e[j];
				}

				//Область Oy_2
				for (int i = 2; i <= i1; i++)
				{
					for (int j = j2 + 1; j <= my - 1; j++) {
						double aUy = abs(Uy[i][j]);
						a[j] = (aUy + Uy[i][j]) / 2. / dy + 1. / Re / dy2;
						c[j] = (aUy - Uy[i][j]) / 2. / dy + 1. / Re / dy2;
						b[j] = 1. / dt + 2. / Re / dy2 + aUy / dy;
						d[j] = vi1[i][j] / dt;
					}

					a[j2] = 0.0;     b[j2] = 1.0;    c[j2] = 0.0;   d[j2] = 2 * (ft[i][j2 + 1] - ft[i][j2]) / dy2;
					a[my] = 0.0;     b[my] = 1.0;    c[my] = 0.0;   d[my] = 2 * (ft[i][my - 1] - ft[i][my]) / dy2;

					Tom(j2, my, a, b, c, d, e);

					for (int j = j2; j <= my; j++)
						vi[i][j] = e[j];
				}

				//Область Oy_3
				for (int i = i1 + 1; i <= mx - 1; i++)
				{
					for (int j = 2; j <= my - 1; j++) {
						double aUy = abs(Uy[i][j]);
						a[j] = (aUy + Uy[i][j]) / 2. / dy + 1. / Re / dy2;
						c[j] = (aUy - Uy[i][j]) / 2. / dy + 1. / Re / dy2;
						b[j] = 1. / dt + 2. / Re / dy2 + aUy / dy;
						d[j] = vi1[i][j] / dt;
					}

					a[1] = 0.0;      b[1] = 1.0;    c[1] = 0.0;      d[1] = 2 * (ft[i][2] - ft[i][1]) / dy2;
					a[my] = 0.0;	 b[my] = 1.0;	c[my] = 0.0;	 d[my] = 2 * (ft[i][my - 1] - ft[i][my]) / dy2;

					Tom(1, my, a, b, c, d, e);

					for (int j = 1; j <= my; j++)
						vi[i][j] = e[j];
				}
			}
			//Границы
			for (int j = 2; j <= my - 1; j++) {
				vi[mx][j] = vi[mx - 1][j];
			}

			for (int j = j1; j <= j2; j++) {
				vi[i1][j] = 2 * (ft[i1 + 1][j] - ft[i1][j]) / dx2;
			}

			for (int i = 2; i <= mx - 1; i++) {
				vi[i][1] = 2 * (ft[i][2] - ft[i][1]) / dy2;
				vi[i][my] = 2 * (ft[i][my - 1] - ft[i][my]) / dy2;
			}

			for (int i = 2; i <= i1; i++) {
				vi[i][j1] = 2 * (ft[i][j1 - 1] - ft[i][j1]) / dy2;
				vi[i][j2] = 2 * (ft[i][j2 + 1] - ft[i][j2]) / dy2;
			}

			for (int i = 1; i <= mx; i++) {
				for (int j = 1; j <= my; j++) {
					vi2[i][j] = vi[i][j] - vis[i][j];
				}
			}
		}

		//Температура
		for(int k = 1; k <= 3; k++ ) //Для ускорения расчетов при больших РR
		{
			//Прогонка по Х			
			{
				//Область Ox_1
				for (int j = 2; j <= j1 - 1; j++) {
					for (int i = 2; i <= mx - 1; i++) {
						double aUx = abs(Ux[i][j]);
						a[i] = (aUx + Ux[i][j]) / 2. / dx + 1. / Re / Pr / dx2;
						c[i] = (aUx - Ux[i][j]) / 2. / dx + 1. / Re / Pr / dx2;
						b[i] = 1. / dt + 2. / Re / Pr / dx2 + aUx / dx;
						d[i] = tets[i][j] / dt;
					}

					a[1] = 0.0;      b[1] = 1.0;    c[1] = 0.0;      d[1] = 0.0;
					a[mx] = 1.0;	 b[mx] = 1.0;	c[mx] = 0.0;	 d[mx] = 0.0;

					Tom(1, mx, a, b, c, d, e);

					for (int i = 1; i <= mx; i++)
						tet1[i][j] = e[i];
				}

				//Область Ox_2
				for (int j = j1; j <= j2; j++) {
					for (int i = i1 + 1; i <= mx - 1; i++) {
						double aUx = abs(Ux[i][j]);
						a[i] = (aUx + Ux[i][j]) / 2. / dx + 1. / Re / Pr / dx2;
						c[i] = (aUx - Ux[i][j]) / 2. / dx + 1. / Re / Pr / dx2;
						b[i] = 1. / dt + 2. / Re / Pr / dx2 + aUx / dx;
						d[i] = tets[i][j] / dt;
					}

					a[i1] = 0.0;  b[i1] = 1.0;       c[i1] = 1.0;       d[i1] = 0.0;
					a[mx] = 1.0;  b[mx] = 1.0;       c[mx] = 0.0;       d[mx] = 0.0;

					Tom(i1, mx, a, b, c, d, e);

					for (int i = i1; i <= mx; i++)
						tet1[i][j] = e[i];
				}

				//Область Ox_3
				for (int j = j2 + 1; j <= my - 1; j++) {
					for (int i = 2; i <= mx - 1; i++) {
						double aUx = abs(Ux[i][j]);
						a[i] = (aUx + Ux[i][j]) / 2. / dx + 1. / Re / Pr / dx2;
						c[i] = (aUx - Ux[i][j]) / 2. / dx + 1. / Re / Pr / dx2;
						b[i] = 1. / dt + 2. / Re / Pr / dx2 + aUx / dx;
						d[i] = tets[i][j] / dt;
					}

					a[1] = 0.0;       b[1] = 1.0;       c[1] = 0.0;     d[1] = 1.0;
					a[mx] = 1.0;	  b[mx] = 1.0;		c[mx] = 0.0;	d[mx] = 0.0;

					Tom(1, mx, a, b, c, d, e);

					for (int i = 1; i <= mx; i++)
						tet1[i][j] = e[i];
				}
			}			

			//Прогонка по Y
			{
				//Область Oy_1
				for (int i = 2; i <= i1; i++)
				{
					for (int j = 2; j <= j1 - 1; j++) {
						double aUy = abs(Uy[i][j]);
						a[j] = (aUy + Uy[i][j]) / 2. / dy + 1. / Re / Pr / dy2;
						c[j] = (aUy - Uy[i][j]) / 2. / dy + 1. / Re / Pr / dy2;
						b[j] = 1. / dt + 2. / Re / Pr / dy2 + aUy / dy;
						d[j] = tet1[i][j] / dt;
					}

					a[1] = 0.0;       b[1] = 1.0;    c[1] = 1.0;        d[1] = 0.0;
					a[j1] = 1.0;	  b[j1] = 1.0;	 c[j1] = 0.0;		d[j1] = 0.0;

					Tom(1, j1, a, b, c, d, e);

					for (int j = 1; j <= j1; j++)
						tet[i][j] = e[j];
				}

				//Область Oy_2
				for (int i = 2; i <= i1; i++)
				{
					for (int j = j2 + 1; j <= my - 1; j++) {
						double aUy = abs(Uy[i][j]);
						a[j] = (aUy + Uy[i][j]) / 2. / dy + 1. / Re / Pr / dy2;
						c[j] = (aUy - Uy[i][j]) / 2. / dy + 1. / Re / Pr / dy2;
						b[j] = 1. / dt + 2. / Re / Pr / dy2 + aUy / dy;
						d[j] = tet1[i][j] / dt;
					}

					a[j2] = 0.0;     b[j2] = 1.0;    c[j2] = 1.0;   d[j2] = 0.0;
					a[my] = 1.0;     b[my] = 1.0;    c[my] = 0.0;   d[my] = 0.0;

					Tom(j2, my, a, b, c, d, e);

					for (int j = j2; j <= my; j++)
						tet[i][j] = e[j];
				}

				//Область Oy_3
				for (int i = i1 + 1; i <= mx - 1; i++)
				{
					for (int j = 2; j <= my - 1; j++) {
						double aUy = abs(Uy[i][j]);
						a[j] = (aUy + Uy[i][j]) / 2. / dy + 1. / Re / Pr / dy2;
						c[j] = (aUy - Uy[i][j]) / 2. / dy + 1. / Re / Pr / dy2;
						b[j] = 1. / dt + 2. / Re / Pr / dy2 + aUy / dy;
						d[j] = tet1[i][j] / dt;
					}

					a[1] = 0.0;     b[1] = 1.0;     c[1] = 1.0;      d[1] = 0.0;
					a[my] = 1.0;	b[my] = 1.0;	c[my] = 0.0;	 d[my] = 0.0;

					Tom(1, my, a, b, c, d, e);

					for (int j = 1; j <= my; j++)
						tet[i][j] = e[j];
				}
			}

			//Границы
			for (int j = 2; j <= my - 1; j++) {
				tet[mx][j] = tet[mx - 1][j];
			}

			for (int i = 2; i <= mx - 1; i++) {
				tet[i][1] = tet[i][2];
				tet[i][my] = tet[i][my - 1];
			}

			for (int i = 2; i <= i1; i++) {
				tet[i][j1] = tet[i][j1 - 1];
				tet[i][j2] = tet[i][j2 + 1];
			}

			for (int j = j1; j <= j2; j++) {
				tet[i1][j] = tet[i1 + 1][j];
			}

			for (int i = 1; i <= mx; i++) {
				for (int j = 1; j <= my; j++) {
					tet2[i][j] = tet[i][j] - tets[i][j];
				}
			}
			
			for (int i = 1; i <= mx; i++) {
				for (int j = 1; j <= my; j++) {
					tets[i][j] = tet[i][j];
				}
			}
		}
		
		for (int i = 1; i <= mx; i++) {
			for (int j = 1; j <= my; j++) {
				vis[i][j] = vi[i][j];
				fts[i][j] = ft[i][j];
				tets[i][j] = tet[i][j];
			}
		}

		if (it % isht == 0) {
			isht += ish;
			sumft = 0.;
			sumvi = 0.;
			sumtet = 0.;

			for (int i = 2; i <= mx - 1; i++) {
				for (int j = 2; j <= my - 1; j++) {
					sumft += abs(ft2[i][j]) / dt;
					sumvi += abs(vi2[i][j]) / dt;
					sumtet += abs(tet2[i][j]) / dt;
				}
			}

			cout << it << "\t" << sumvi << "\t" << sumtet << endl;
		}

		if (it < itk) continue;

		cout << "Ввод дополнительных итераций" << endl;
		int itdop;
		cin >> itdop;
		itk += itdop;

		if (itdop > 0) {
			isht -= ish;
			cout << "Ввод шага сходимости" << endl;
			cin >> ish;
			isht += ish;
		}

		else if (itdop == 0) break;
	}




	//Сохранение в файл
	cout << "Введите имя варианта решения" << endl;
	string name;
	cin >> name;

	//Скорости
	ofstream fout_izo_tem_ft(name + "_izo_tem_ft.dat");
	for (int i = 1; i <= mx; i++) {
		for (int j = 1; j <= my; j++) {
			x = (i - 1) * dx;
			y = (j - 1) * dy;
			fout_izo_tem_ft << x <<" " << y << " " << tet[i][j] << " " << ft[i][j] << endl;
		}
	}
	fout_izo_tem_ft.close();
	
	ofstream fout_izo_ux_uy_vi(name + "_izo_ux_uy_vi.dat");
	for (int i = 1; i <= mx; i++) {
		for (int j = 1; j <= my - 1; j++) {
			x = (i - 1) * dx;
			y = (j - 1) * dy;
			fout_izo_ux_uy_vi << x <<" " << y << " " << Ux[i][j] << " " << Uy[i][j] << " " << vi[i][j] << endl;
		};
	}
	fout_izo_ux_uy_vi.close();

	//Ux по сечениям
	ofstream fout_Ux_nizhnee_sechenie(name + "_Ux_nizhnee_sechenie.dat");	
		for (int j = 1; j <= j1; j++) {
			y = (j - 1) * dy;
			fout_Ux_nizhnee_sechenie << y << " " << Ux[2][j] << " " << Ux[i1/4][j] << " " << Ux[i1/2][j] << " " << Ux[i1-1][j] << endl;
		}	
	fout_Ux_nizhnee_sechenie.close();

	ofstream fout_Ux_verhnee_sechenie(name + "_Ux_verhnee_sechenie.dat");	
		for (int j = j2; j <= my; j++) {
			y = (j - 1) * dy;
			fout_Ux_verhnee_sechenie << y << " " << Ux[2][j] << " " << Ux[i1/4][j] << " " << Ux[i1/2][j] << " " << Ux[i1-1][j] << endl;
		}
		fout_Ux_verhnee_sechenie.close();
	
	ofstream fout_Ux_obshaia_chast(name + "_Ux_obshaia_chast.dat");	
		for (int j = 1; j <= my; j++) {
			y = (j - 1) * dy;
			fout_Ux_obshaia_chast << y << " " << Ux[(i1+mx)/2][j] << " " << Ux[2 * (i1 + mx) / 3][j] << " " << Ux[mx][j] << endl;
		}
	fout_Ux_obshaia_chast.close();
	
	//Температура по сечениям
	ofstream fout_tet_nizhnee_sechenie(name + "_tet_nizhnee_sechenie.dat");
	for (int j = 1; j <= j1; j++) {
		y = (j - 1) * dy;
		fout_tet_nizhnee_sechenie << y << " "  << tet[i1 / 2][j] << " " << tet[i1 - 1][j] << endl;
	}
	fout_Ux_nizhnee_sechenie.close();

	ofstream fout_tet_verhnee_sechenie(name + "_tet_verhnee_sechenie.dat");
	for (int j = j2; j <= my; j++) {
		y = (j - 1) * dy;
		fout_tet_verhnee_sechenie << y << " " << tet[i1 / 2][j] << " " << tet[i1 - 1][j] << endl;
	}
	fout_tet_verhnee_sechenie.close();

	ofstream fout_tet_obshaia_oblast(name + "_tet_obshaia_oblast.dat");
	for (int j = 1; j <= my; j++) {
		y = (j - 1) * dy;
		fout_tet_obshaia_oblast << y << " " << tet[(i1 + mx) / 2][j] << " " << tet[2 * (i1 + mx) / 3][j] << " " << tet[mx][j] << endl;
	}
	fout_tet_obshaia_oblast.close();

	//Uy по сечениям
	ofstream fout_Uy_nizhnee_sechenie(name + "_Uy_nizhnee_sechenie.dat");
	for (int j = 1; j <= j1; j++) {
		y = (j - 1) * dy;
		fout_Uy_nizhnee_sechenie << y << " " << Ux[i1 / 2][j] << " " << Ux[i1 - 1][j] << endl;
	}
	fout_Uy_nizhnee_sechenie.close();

	ofstream fout_Uy_verhnee_sechenie(name + "_Uy_verhnee_sechenie.dat");
	for (int j = j2; j <= my; j++) {
		y = (j - 1) * dy;
		fout_Uy_verhnee_sechenie << y << " " << Ux[i1 / 2][j] << " " << Ux[i1 - 1][j] << endl;
	}
	fout_Uy_verhnee_sechenie.close();

	ofstream fout_Uy_obshaia_oblast(name + "_Uy_obshaia_oblast.dat");
	for (int j = 1; j <= my; j++) {
		y = (j - 1) * dy;
		fout_Uy_obshaia_oblast << y << " " << Uy[(i1 + mx) / 2][j] << " " << Uy[2 * (i1 + mx) / 3][j] << " " << Uy[mx][j] << endl;
	}
	fout_Uy_obshaia_oblast.close();
}

