#include <stdio.h>
#include <math.h>

#define IMM 35
#define JMM 21
#define PI atan(1)*4
const int IM = IMM;
const int JM = JMM;
const double CHORD = 1.0;
const double RADIUS = 1.5;
const double T = 0.2;
const double ERRMAX = 0.01;

double data[2][JMM][IMM];

double a[JMM - 2][IMM];
double b[JMM - 2][IMM];
double c[JMM - 2][IMM];

double above[IMM - 3];
double below[IMM - 3];
double diagonal[IMM - 2];
double constant[IMM - 2];
double solution[IMM - 2];

double airfoil(double x)
{
    return (T / 0.2) * ((0.2969 * sqrt(x)) - (0.126 * x) - (0.3516 * pow(x, 2.0)) + (0.2843 * pow(x, 3.0)) - (0.1015 * pow(x, 4.0)));
}

void init()
{
    const int NOSE = ((IM + 1) / 2) - 1;
    const double DELTAX = CHORD / NOSE;
    const double DELTAI = (2.0 * PI) / (IM - 1);


    data[0][0][0] = CHORD / 2;
    data[1][0][0] = 0;

    data[0][0][NOSE] = -CHORD / 2;
    data[1][0][NOSE] = 0;

    data[0][0][IM - 1] = CHORD / 2;
    data[1][0][IM - 1] = 0;

    for (int i = 1; i < NOSE; i++)
    {
        data[0][0][NOSE + i] = (i * DELTAX * CHORD) - (CHORD / 2);
        data[0][0][NOSE - i] = data[0][0][NOSE + i];
        data[1][0][NOSE + i] = airfoil(i * DELTAX * CHORD);
        data[1][0][NOSE - i] = (-1.0) * data[1][0][NOSE + i];
    }

    for (int i = 0; i < IM; i++)
    {
        data[0][JM - 1][i] = RADIUS * cos(DELTAI * i);
        data[1][JM - 1][i] = -RADIUS * sin(DELTAI * i);
    }

    double delx;
    double dely;

    for (int i = 0; i < IM; i++)
    {
        delx = (data[0][JM - 1][i] - data[0][0][i]) / (JM - 1);
        dely = (data[1][JM - 1][i] - data[1][0][i]) / (JM - 1);

        for (int j = 1; j < JM - 1; j++)
        {
            data[0][j][i] = (double)j * delx + data[0][0][i];
            data[1][j][i] = (double)j * dely + data[1][0][i];
        }
    }
}

void initCoef()
{
    double xEta;
    double xXi;
    double yEta;
    double yXi;

    for (int j = 1; j < JM - 1; j++)
    {

        xEta = (data[0][j + 1][0] - data[0][j - 1][0]) / 2;
        yEta = (data[1][j + 1][0] - data[1][j - 1][0]) / 2;
        xXi = (data[0][j][1] - data[0][j][IM - 2]) / 2;
        yXi = (data[1][j][1] - data[1][j][IM - 2]) / 2;

        a[j - 1][0] = (xEta * xEta) + (yEta * yEta);
        b[j - 1][0] = (xEta * xXi) + (yEta * yXi);
        c[j - 1][0] = (xXi * xXi) + (yXi * yXi);

    }

    for (int i = 1; i < IM - 1; i++)
    {
        for (int j = 1; j < JM - 1; j++)
        {
            xEta = (data[0][j + 1][i] - data[0][j - 1][i]) / 2;
            yEta = (data[1][j + 1][i] - data[1][j - 1][i]) / 2;
            xXi = (data[0][j][i + 1] - data[0][j][i - 1]) / 2;
            yXi = (data[1][j][i + 1] - data[1][j][i - 1]) / 2;

            a[j - 1][i] = (xEta * xEta) + (yEta * yEta);
            b[j - 1][i] = (xEta * xXi) + (yEta * yXi);
            c[j - 1][i] = (xXi * xXi) + (yXi * yXi);
        }
    }

    for (int j = 1; j < JM - 1; j++)
    {

        xEta = (data[0][j + 1][IM - 1] - data[0][j - 1][IM - 1]) / 2;
        yEta = (data[1][j + 1][IM - 1] - data[1][j - 1][IM - 1]) / 2;
        xXi = (data[0][j][1] - data[0][j][IM - 2]) / 2;
        yXi = (data[1][j][1] - data[1][j][IM - 2]) / 2;

        a[j - 1][IM - 1] = (xEta * xEta) + (yEta * yEta);
        b[j - 1][IM - 1] = (xEta * xXi) + (yEta * yXi);
        c[j - 1][IM - 1] = (xXi * xXi) + (yXi * yXi);

    }
}

void initThomasX(int j)
{
    for (int i = 1; i < IM - 2; i++)
    {
        above[i - 1] = a[j - 1][i];
        below[i - 1] = a[j - 1][i + 1];
    }

    for (int i = 1; i < IM - 1; i++)
    {
        diagonal[i - 1] = (-2.0) * (a[j - 1][i] + c[j - 1][i]);
        constant[i - 1] = ((b[j - 1][i] / 2.0) * (data[0][j + 1][i + 1] - data[0][j - 1][i + 1] + data[0][j - 1][i - 1] - data[0][j + 1][i - 1])) - (c[j - 1][i] * (data[0][j + 1][i] + data[0][j - 1][i]));
    }

    constant[0] = constant[0] - (a[j - 1][1] * data[0][j][0]);
    constant[IM - 3] = constant[IM - 3] - (a[j - 1][IM - 2] * data[0][j][IM - 1]);
}

void initThomasY(int j)
{
    for (int i = 1; i < IM - 2; i++)
    {
        above[i - 1] = a[j - 1][i];
        below[i - 1] = a[j - 1][i + 1];
    }

    for (int i = 1; i < IM - 1; i++)
    {
        diagonal[i - 1] = (-2.0) * (a[j - 1][i] + c[j - 1][i]);
        constant[i - 1] = ((b[j - 1][i] / 2.0) * (data[1][j + 1][i + 1] - data[1][j - 1][i + 1] + data[1][j - 1][i - 1] - data[1][j + 1][i - 1])) - (c[j - 1][i] * (data[1][j + 1][i] + data[1][j - 1][i]));
    }

    constant[0] = constant[0] - (a[j - 1][1] * data[1][j][0]);
    constant[IM - 3] = constant[IM - 3] - (a[j - 1][IM - 2] * data[1][j][IM - 1]);
}

void thomas()
{
    for (int i = 1; i < IM - 2; i++)
    {
        diagonal[i] = diagonal[i] - ((below[i - 1] * above[i - 1]) / diagonal[i - 1]);
        constant[i] = constant[i] - ((constant[i - 1] * below[i - 1]) / diagonal[i - 1]);
    }

    solution[IM - 3] = constant[IM - 3] / diagonal[IM - 3];

    for (int i = IM - 4; i > -1; i--)
    {
        solution[i] = (constant[i] - (above[i] * solution[i + 1])) / diagonal[i];
    }
}

double adjustX()
{
    double temp = 0;

    for (int j = 1; j < JM - 1; j++)
    {
        data[0][j][0] = a[j - 1][0] * (data[0][j][1] + data[0][j][IM - 2]);
        data[0][j][0] = data[0][j][0] + c[j - 1][0] * (data[0][j + 1][0] + data[0][j - 1][0]);
        data[0][j][0] = data[0][j][0] - (b[j - 1][0] / 2) * (data[0][j + 1][1] - data[0][j - 1][1] + data[0][j + 1][IM - 2] - data[0][j - 1][IM - 2]);
        data[0][j][0] = data[0][j][0] / (2 * (a[j - 1][0] + c[j - 1][0]));

        temp = temp + fabs(data[0][j][0] - data[0][j][IM - 1]);

        data[0][j][IM - 1] = data[0][j][0];
    }

    return temp;
}

double adjustY()
{
    double temp = 0;

    for (int j = 1; j < JM - 1; j++)
    {
        data[1][j][0] = a[j - 1][0] * (data[1][j][1] + data[1][j][IM - 2]);
        data[1][j][0] = data[1][j][0] + c[j - 1][0] * (data[1][j + 1][0] + data[1][j - 1][0]);
        data[1][j][0] = data[1][j][0] - (b[j - 1][0] / 2) * (data[1][j + 1][1] - data[1][j - 1][1] + data[1][j + 1][IM - 2] - data[1][j - 1][IM - 2]);
        data[1][j][0] = data[1][j][0] / (2 * (a[j - 1][0] + c[j - 1][0]));

        temp = temp + fabs(data[1][j][0] - data[1][j][IM - 1]);

        data[1][j][IM - 1] = data[1][j][0];
    }

    return temp;
}

void print(const char* name)
{
    FILE* file;
    file = fopen(name, "w");

    for (int j = 0; j < JM; j++)
    {
        for (int i = 0; i < IM; i++)
        {
            fprintf(file, "%.6f %.6f\n", data[0][j][i], data[1][j][i]);
        }
    }

    fclose(file);
}

int main()
{

    double ERRX;
    double ERRY;

    init();

    print("AlgerbraicGrid.dat");

    do
    {
        ERRX = 0;
        ERRY = 0;

        initCoef();

        for (int j = 1; j < JM - 1; j++)
        {
            initThomasX(j);
            thomas();

            for (int i = 1; i < IM - 1; i++)
            {
                ERRX = ERRX + fabs(solution[i - 1] - data[0][j][i]);
                data[0][j][i] = solution[i - 1];
            }
            ERRX = ERRX + adjustX();

            initThomasY(j);
            thomas();

            for (int i = 1; i < IM - 1; i++)
            {
                ERRY = ERRY + fabs(solution[i - 1] - data[1][j][i]);
                data[1][j][i] = solution[i - 1];
            }
            ERRY = ERRY + adjustY();

        }

    } while ((ERRX + ERRY) > ERRMAX);

    print("EllipticGrid.dat");

    return 0;
}