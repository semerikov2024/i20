// реалізація лінійної регресії

using namespace std; // для використання стандартних бібліотек

#include <iostream> // для роботи з потоками
#include <fstream> // для роботи з файловими потоками
#include <vector> // для роботи з векторами
#include <string> // для роботи з рядками
#include <math.h> // для роботи з sqrt, ...

// дійсний вектор
typedef vector<double> dvec;
// дійсна матриця
typedef vector<dvec> dmat;

// оператор введення вектору з потоку
istream& operator>>(istream &s, dvec &v)
{
    int size = v.size(); // розмірність вектору

    for(int i=0;i<size;i++) // ітерація
        s>>v[i]; // читання одного елемента з потоку
    return s;
}


// оператор виведення вектору з потоку
ostream& operator<<(ostream &s, dvec v)
{
    int size = v.size(); // розмірність вектору

    for(int i=0;i<size;i++) // ітерація
        s<<v[i]<<"\t"; // запис одного елемента у потік
    return s;
}

// дії над векторами

dvec operator+(dvec a, dvec b)// додавання: c=a+b
{
    // перевірка розмірностей 
    if(a.size() != b.size())
        throw string("Розмірності мають бути однаковими");
    dvec result(a.size());

    for(int i=0;i<result.size();i++)
        result[i]=a[i]+b[i];

    return result;
}

dvec operator+(dvec a)// збереження знаку: c = +a
{
    return a;
}


dvec operator-(dvec a) // зміна знаку: c = -a
{
    dvec result(a.size());

    for(int i=0;i<result.size();i++)
        result[i]=-a[i];

    return result;
}


dvec operator-(dvec a, dvec b) // оператор віднімання c = a - b
{
    // a - b == a + (-b)
    return a + (-b);
}

dvec operator*(dvec a, double n)// множення вектора на число: c=a*n
{
    dvec result(a.size());

    for(int i=0;i<result.size();i++)
        result[i]=a[i]*n;

    return result;
}

dvec operator*(double n, dvec a)// множення числа на вектор : n * a == a * n
{
    return a * n;
}

double operator*(dvec a, dvec b)// скалярний добуток: c=a*b
{
    // перевірка розмірностей 
    if(a.size() != b.size())
        throw string("Розмірності мають бути однаковими");
    double result=0;

    for(int i=0;i<a.size();i++)
        result+=a[i]*b[i];

    return result;
}

double abs(dvec a) // vector length = sqrt(a*a)
{
    return sqrt(a*a);
}



// матриця
class matrix
{
    private:
    // дані
    int rows, // кількість рядків
        cols; // кількість стовпців
    dmat mtr; // масив для зберігання даних

    // методи
    public:
    // дані

    // методи
    matrix(int m, int n, double filler); // конструктор матриці певного розміру
    matrix(const matrix&); //конструктор копіювання
    dvec& operator[](int);// індексація матриці [] - рядок
    dvec operator()(int);// індексація матриці [] - стовпець
    matrix operator=(matrix); // оператор надання значення =
    matrix operator~(); // оператор транспонування ~
    //
    int getrows() {return rows;}
    int getcols() {return cols;}
    //...
    friend istream& operator>>(istream &s, matrix &m);
    friend ostream& operator<<(ostream &s, matrix m);
};


// конструктор матриці певного розміру
matrix::matrix(int m=1, int n=1, double filler=0): rows(m), cols(n)
{
    // перевірка розмірностей 
    if(m<1 || n<1)
        throw string("Розмірності мають бути натуральними числами");
    mtr=dmat(rows); // матриця заданої розмірності у рядках з порожніми рядками

    //розширення кожного рядка до заданої кількості стовпців
    for(int i=0;i<rows;i++)
    {
        mtr[i]=dvec(cols);
        for(int j=0;j<cols;j++)
            mtr[i][j] = filler;
    }
}

//конструктор копіювання
matrix::matrix(const matrix &exist): rows(exist.rows), cols(exist.cols)
{
    mtr=dmat(rows); // матриця заданої розмірності у рядках з порожніми рядками

    //розширення кожного рядка до заданої кількості стовпців
    for(int i=0;i<rows;i++)
    {
        mtr[i]=dvec(cols);
        for(int j=0;j<cols;j++)
            mtr[i][j] = exist.mtr[i][j];
    }
}


//оператор надання значення
matrix matrix::operator=(matrix exist)
{
    rows=exist.rows, cols=exist.cols;

    mtr=dmat(rows); // матриця заданої розмірності у рядках з порожніми рядками

    //розширення кожного рядка до заданої кількості стовпців
    for(int i=0;i<rows;i++)
    {
        mtr[i]=dvec(cols);
        for(int j=0;j<cols;j++)
            mtr[i][j] = exist.mtr[i][j];
    }

    return *this;
}


dvec& matrix::operator[](int row)
{
    // перевірка розмірностей 
    if(row < 0 || row >= rows)
    //if(this->row < 0 || (*this).row >= rows)
        throw string("Номер рядка має бути натуральним числом та не перевищувати кількості рядків");
    return mtr[row];
    //return this->mtr[row];
}

dvec matrix::operator()(int col) // індексація матриці () - стовпець
{
    // перевірка розмірностей 
    if(col < 0 || col >= cols)
        throw string("Номер стовпця має бути натуральним числом та не перевищувати кількості стовпців");
    dvec result(rows);
    
    for(int i=0;i<rows;i++)
        result[i]=mtr[i][col];
    return result;
}


// оператор введення матриці з потоку
istream& operator>>(istream &s, matrix &m)
{
    int size = m.rows; // розмірність матриці == кількість рядків

    for(int i=0;i<size;i++) // ітерація
        s>>m.mtr[i]; // читання одного елемента з потоку
    return s;
}


// оператор виведення матриці у потік
ostream& operator<<(ostream &s, matrix m)
{
    int size = m.rows; // розмірність матриці == кількість рядків

    for(int i=0;i<size;i++) // ітерація
        s<<m.mtr[i]<<"\n"; // запис одного елемента у потік
    return s;
}


// оператор транспонування ~
matrix matrix::operator~()
{
    matrix T(cols, rows);

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            T[j][i] = mtr[i][j];
    return T;
}


// дії над матрицями

matrix operator+(matrix a, matrix b)// додавання: c=a+b
{
    // перевірка розмірностей 
    if(a.getrows() != b.getrows() || a.getcols() != b.getcols())
        throw string("Розмірності мають бути однаковими");
    matrix result(a.getrows(), a.getcols());

    for(int i=0;i<result.getrows();i++)
        result[i]=a[i]+b[i];

    return result;
}

matrix operator+(matrix a)// збереження знаку: c = +a
{
    return a;
}


matrix operator-(matrix a) // зміна знаку: c = -a
{
    matrix result = a;

    for(int i=0;i<result.getrows();i++)
        result[i]=-a[i];

    return result;
}


matrix operator-(matrix a, matrix b) // оператор віднімання c = a - b
{
    // a - b == a + (-b)
    return a + (-b);
}

matrix operator*(matrix a, double n)// множення вектора на число: c=a*n
{
    matrix result = a;

    for(int i=0;i<result.getrows();i++)
        result[i]=a[i]*n;

    return result;
}

matrix operator*(double n, matrix a)// множення числа на вектор : n * a == a * n
{
    return a * n;
}

matrix operator*(matrix a, matrix b)// добуток: c=a*b
{
    /*
    a=[r1, c1] 
    1 2 3
    4 5 6

    b=[r2, c2]
    7 8 3 6
    9 0 4 7
    1 2 5 -1

    c=[r1, c2]
    a[1]*b(1) a[1]*b[S2] a[1]*b[S3] a[1]*b[S4] 
    a[2]*b[S1] a[2]*b[S2] a[2]*b[S3] a[2]*b[S4] 
    */

   int r1 = a.getrows(), c1 = a.getcols(),
        r2 = b.getrows(), c2 = b.getcols();

    // перевірка розмірностей 
    if(c1 != r2)
        throw string("Множення матриць неможливо через несумірність");
    matrix result(r1, c2);

    for(int i=0;i<r1;i++)
        for(int j=0;j<c2;j++)
            result[i][j]=a[i]*b(j);

    return result;
}


matrix operator^(matrix a, int n);

// знаходження оберненої матриці
matrix inverse(matrix a)
{
    // перевірка розмірностей 
    if(a.getrows() != a.getcols())
        throw string("Матриця не квадратна");

    // M = [A | E] - побудова розширеної матриці
    matrix E = a^0;
    matrix M(a.getrows(), a.getcols()*2);
    for(int i=0;i<a.getrows();i++)
        for(int j=0;j<a.getcols();j++)
        {
            M[i][j]=a[i][j];
            M[i][j+a.getcols()]=E[i][j];
        }

    //прямий хід Гауса (нулі під головною діагоналлю)
    for(int i=0;i<M.getrows();i++)
        for(int j=i+1;j<M.getrows();j++)
        {
            double c = - M[j][i] / M[i][i]; // головний елемент
            // Додавання до одного рядка іншого рядка помноженого на скаляр.
            M[j] = M[j] + M[i] * c;
        }

    //обернений хід Гауса (нулі над головною діагоналлю)

    for(int i=M.getrows()-1;i>=0;i--)
        for(int j=i-1;j>=0;j--)
        {
            double c = - M[j][i] / M[i][i]; // головний елемент
            // Додавання до одного рядка іншого рядка помноженого на скаляр.
            M[j] = M[j] + M[i] * c;
        }

    //одиниці на головній діагоналі
    for(int i=0;i<M.getrows();i++)
    {
        double c = 1/M[i][i];
        // Множення рядка на не нульове скалярне значення.
        M[i] = M[i] * c;
    }

    matrix res=a;
    // [E | A^-1]
    for(int i=0;i<a.getrows();i++)
        for(int j=0;j<a.getcols();j++)
            res[i][j] = M[i][j+a.getcols()];

    return res;
}



matrix operator^(matrix a, int n) // a^n
{
    // перевірка розмірностей 
    if(a.getrows() != a.getcols())
        throw string("Матриця не квадратна");

    if(n==1) 
        return a;
    else
        if(n==0)
        {
            matrix result(a.getrows(), a.getcols());
            for(int i=0;i<result.getrows();i++)
                result[i][i] = 1;
            return result;
        }
        if(n>1)
            return a * a^(n-1); // a^n == a * a^(n-1)
        else
            if(n==-1)
                return inverse(a);
            else
                // a^(-n) == (a^(-1))^n
                return inverse(a)^(-n);
}


double abs(matrix a) // vector length = sqrt(a*a)
{
    return 0;//sqrt(a*a);
}



// формування матриць X, Y
void getXY(matrix &data, matrix &X, int numx, matrix &Y, int numy)
{
    int rows=data.getrows();

    X=matrix(rows, numx+1, 1); // x1	x2	x3	x4	x5	x6	x7	x8	1
    Y=matrix(rows, numy); // y - вік

    for(int i=0; i<rows; i++)
    {
        for(int j=0;j<numx;j++)
            X[i][j] = data[i][j];

        int startcolY=data.getcols()-numy;
        for(int j=0;j<numy;j++)
            Y[i][j] = data[i][startcolY+j]+1.5;
    }
}


int main()
{
    ifstream f("data/abalone_tabular.txt"); // спроба відкриття файлу
    if(!f) 
    {
        cerr<<"Неможливо відкрити файл data/abalone_tabular.txt"<<endl;
        return 1;
    }

    int num_Instances,  // кількість вимірювань (рядків)
        num_Features_x, // кількість незалежних змінних
        num_Features_y; // кількість залежних змінних

    // читання розмірностей з файлу
    if(!( ((f>>num_Instances)>>num_Features_x)>>num_Features_y ))
    {
        cerr<<"Помилка читання з файлу data/abalone_tabular.txt"<<endl;
        return 2;
    }

    // перевірка розмірностей на коректність
    if( (num_Instances < 1) || (num_Features_x < 1) || (num_Features_y < 1) )
    {
        cerr<<"Розмірності мають бути натуральними числами"<<endl;
        return 3;
    }

    cout<<"Kількість вимірювань (рядків) " << num_Instances
        <<", кількість незалежних змінних " << num_Features_x
        <<", кількість залежних змінних " << num_Features_y
        <<endl;

    matrix data(num_Instances, num_Features_x+num_Features_y);
    f>>data;
    matrix X, Y;
    getXY(data, X, num_Features_x, Y, num_Features_y);

    //cout<<"X"<<endl<<X<<endl;
    //cout<<"Y"<<endl<<Y<<endl;

    matrix XT = ~X; // транспонування
    //matrix C = XT * X;
    //cout<<"XT*X"<<endl<<C<<endl;
    matrix A=((XT*X)^(-1))*XT*Y;
    cout<<"A"<<endl<<A<<endl;

/*
    dvec a=X[2];
    dvec b=X[0];
    cout<<"a   = "<<a<<endl;
    cout<<"b   = "<<b<<endl;
    dvec c=a-b;
    cout<<"a-b = "<<c<<endl;
    cout<<"c*2 = "<<c*2<<endl;
    cout<<"a*b = "<<a*b<<endl;
    cout<<"|a| = "<<abs(a)<<endl;
*/

    //cout<<"Воно працює"<<endl;
    return 0;
}
