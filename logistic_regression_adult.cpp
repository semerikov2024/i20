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
            Y[i][j] = (data[i][startcolY+j]==2) ? 1 : 0;
    }
}


// формування матриць X [1, 2 pow], Y
void getX2Y(matrix &data, matrix &X, int numx, matrix &Y, int numy)
{
    int rows=data.getrows();

    X=matrix(rows, numx*2+1, 1); // x1	x2	x3	x4	x5	x6	x7	x8	1
    Y=matrix(rows, numy); // y - вік

    for(int i=0; i<rows; i++)
    {
        for(int j=0;j<numx;j++)
        {
            X[i][j*2] = data[i][j];
            X[i][j*2+1] = pow(data[i][j],2);
        }

        int startcolY=data.getcols()-numy;
        for(int j=0;j<numy;j++)
            Y[i][j] = (data[i][startcolY+j]==2) ? 1 : 0;
    }
}


// logistic
double sigma(double x, double a=1, double b=1, double c=1, double d=0)
{
    return c/(1+b*exp(-a*x))+d;
}

// визначення y^
double func_y_hat(matrix A, dvec x)
{
    return sigma(A(0)*x);
}

// визначення J (помилка = logistic regression)
double J(matrix A, matrix X, matrix Y)
{
    double result = 0;
    int m = X.getrows();

    for(int i=0;i<m;i++)
        result += Y[i][0] * log(func_y_hat(A, X[i])) + (1 - Y[i][0]) * log(1 - func_y_hat(A, X[i]));

    return -(1.0/m)*result;
}

// градієнтний спуск
matrix grad_descent(matrix X, matrix Y, string str, double alpha=0.01, int epochs = 10, double epsilon=0.0001, bool out=false)
{
    matrix A(X.getcols(), 1); // w1, w2, ... , wn, b;
    int n=X.getcols(),
        m=X.getrows();

    matrix grad = A;
    double error_start, error_end;

    do
    {
        error_start = J(A, X, Y);
        for(int k=0;k<n;k++)
        {
            grad[k][0]=0;
            for(int i=0;i<m;i++)
                grad[k][0] += (func_y_hat(A, X[i]) - Y[i][0]) * X[i][k];
        }
        A = A - (alpha/m) * grad;
        error_end = J(A, X, Y);
        if(error_end<error_start)
            break;
        else
        {
            A = matrix(X.getcols(), 1);
            alpha *= 0.9;
            if(out)
                cout<<"Крок скориговано у "<<alpha<<endl;
        }
    }
    while(error_end>error_start);

    ifstream fweights(str);
    if(fweights)
    {
        int A_g_rows, A_g_cols;
        fweights>>A_g_rows>>A_g_cols;
        fweights>>A;
    }

    for(int num=0;num<epochs;num++)
    {
        error_start = error_end;
        for(int k=0;k<n;k++)
        {
            grad[k][0]=0;
            for(int i=0;i<m;i++)
                grad[k][0] += (func_y_hat(A, X[i]) - Y[i][0]) * X[i][k];
        }
        A = A - (alpha/m) * grad;
        error_end = J(A, X, Y);
        if(out)
            cout<<"epoch = "<<(num+1)<<", J = " << error_end <<endl;
        if(abs(error_end - error_start) < epsilon)
            break;
        if(error_end>error_start)
            alpha *= 0.9;
        if(num%100==0)
        {
            ofstream file(str);
            file<<A.getrows()<<" "<<A.getcols()<<endl<<A<<endl;
        }
    }

    ofstream file(str);
    file<<A.getrows()<<" "<<A.getcols()<<endl<<A<<endl;

    return A;
}

// normalize by mean
matrix norm_mean(matrix X)
{
    int m = X.getrows();
    matrix X_n=X;

    for(int i=0;i<X.getcols()-1;i++)
    {
        double min_col=X[0][i], max_col=X[0][i], ave_col=0;
        for(int j=0;j<X.getrows();j++)
        {
            ave_col += X[j][i];
            if(min_col>X[j][i])
                min_col = X[j][i];
            if(max_col<X[j][i])
                max_col = X[j][i];
        }
        ave_col /= m;
        for(int j=0;j<X_n.getrows();j++)
            X_n[j][i] = (X[j][i] - ave_col) / (max_col - min_col);
    }
    return X_n;
}


int main()
{
    ifstream f("data/breast-cancer-wisconsin_tabular.txt"); // спроба відкриття файлу
    if(!f) 
    {
        cerr<<"Неможливо відкрити файл data/breast-cancer-wisconsin_tabular.txt"<<endl;
        return 1;
    }

    int num_Instances,  // кількість вимірювань (рядків)
        num_Features_x, // кількість незалежних змінних
        num_Features_y; // кількість залежних змінних

    // читання розмірностей з файлу
    if(!( ((f>>num_Instances)>>num_Features_x)>>num_Features_y ))
    {
        cerr<<"Помилка читання з файлу data/breast-cancer-wisconsin_tabular.txt"<<endl;
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

    //getXY(data, X, num_Features_x, Y, num_Features_y);
    getX2Y(data, X, num_Features_x, Y, num_Features_y);

  //  cout<<"X:"<<endl<<X<<endl;

    matrix X_n = norm_mean(X);

   // cout<<"X_n:"<<endl<<X_n<<endl;

/*
    matrix XT = ~X_n; // транспонування
    matrix A=((XT*X_n)^(-1))*XT*Y;
    cout<<"A (точне)"<<endl<<A<<endl;
    cout<<"J = "<<J(A, X_n, Y)<<endl;
*/

    matrix A_g=grad_descent(X_n, Y, "weights2.txt", 1, 100, 0.000001, true);

/*
    ifstream fweights("weights.txt");
    int A_g_rows, A_g_cols;
    fweights>>A_g_rows>>A_g_cols;
    matrix A_g(A_g_rows, A_g_cols);
    fweights>>A_g;
    cout<<"A (grad)"<<endl<<A_g<<endl;
    cout<<"J = "<<J(A_g, X, Y)<<endl;
*/

    int count = 0, count_1 = 0, count_0 = 0, count_0_instead_1 = 0, count_1_instead_0 = 0;
    for(int i=0;i<Y.getrows();i++)
    {
        double Y_hat = func_y_hat(A_g, X_n[i]); 
        
        if(Y[i][0])
            count_1++;
        else 
            count_0++;
        double sol=(Y_hat >= 0.5) ? 1 : 0;
        if(Y[i][0] && !sol)
            count_0_instead_1++;
        if(!Y[i][0] && sol)
            count_1_instead_0++;
        cout<<"Для i = "<<(i+1)<<" y = "<<Y[i][0]<<", "<<
            "Y_hat = " << Y_hat << ", P = " << sol <<
            ((Y[i][0] == sol) ? ", OK" : ", ERROR")
            << endl;
        if(Y[i][0] != sol)
            count++;
    }    
    cout<<endl<<"Error % = "<< (100.*count/Y.getrows())<<endl;
    cout<<endl<<"Error 1 instead 0 (false positive) % = "<< (100.*count_1_instead_0/count_0)<<endl;
    cout<<endl<<"Error 0 instead 1 (false negative) % = "<< (100.*count_0_instead_1/count_1)<<endl;

    //cout<<"Воно працює"<<endl;
    return 0;
}
