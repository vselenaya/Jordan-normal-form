// =================================================================================================
// АЛГОРИТМ ПОИСКА ЖОРДАНОВОЙ НОРМАЛЬНОЙ ФОРМЫ (ЖНФ) МАТРИЦЫ
// =================================================================================================
//
// ЦЕЛЬ:
// Для заданной квадратной матрицы A (размера n x n) и её вещественных собственных чисел
// (с учётом их алгебраических кратностей, те кратностей как корней характеристического многочлена) найти:
// 1. Жорданову нормальную форму J – блочно-диагональную матрицу, где каждый блок
//    (жорданова клетка) соответствует одному собственному значению.
// 2. Матрицу перехода S (матрицу смены базиса) такую, что A = S * J * S^(-1).
//    Столбцы матрицы S образуют базис (жорданов базис), в котором оператор,
//    задаваемый матрицей A, имеет вид J.
//
// В данном случае мы считаем, что у матрицы есть все n собственных значений, поэтому Жорданова форма существует. В общем
// случае, если это не так, можно расширить поле (над которым ищем корни характеристического многочлена матрицы) -- например,
// поле вещественных чисел можно расширить до комплексного.
//
//
// ОСНОВНАЯ ИДЕЯ АЛГОРИТМА:
//
// Алгоритм строится на основе теории корневых подпространств и жордановых цепочек.
//
// 1. КОРНЕВЫЕ ПОДПРОСТРАНСТВА:
//    Для каждого собственного значения λ (с алгебраической кратностью k_λ) существует
//    корневое подпространство K_λ := Ker((A - λI)^m), где m – такое наименьшее натуральное
//    число, что пространство стабилизируется: Ker((A - λI)^m) = Ker((A - λI)^(m+1)) = Ker((A - λI)^(m+2)) = ...
//    [это то же самое, что размерность совпала с кратностью: dim(Ker((A - λI)^m)) = k_λ].
//    Такое число m называется высотой корневого вектора.
//    Всё пространство V (в нашем случае R^n, где R - поле вещественных чисел) является прямой суммой корневых
//    подпространств: V = K_λ1 ⊕ K_λ2 ⊕ ... ⊕ K_λs.
//    Задача сводится к построению жорданова базиса для каждого корневого
//    подпространства K_λ независимо. Тогда матрица оператора A будет блочно-диагональной,
//    причем i-й диагональный блок будет матрицей ограничения A на подпространство K_λi. Это и есть Жорданова
//    нормальная форма J.
//
// 2. ЖОРДАНОВЫ ЦЕПОЧКИ:
//    Внутри каждого корневого подпространства K_λ базис строится из жордановых цепочек.
//    Жорданова цепочка длины p, порождённая вектором v_p (называемым присоединённым
//    вектором порядка p или "вершиной" цепочки), имеет вид:
//      v_1 = (A - λI)^(p-1) * v_p   (собственный вектор)
//      v_2 = (A - λI)^(p-2) * v_p
//      ...
//      v_{p-1} = (A - λI) * v_p
//      v_p                            (присоединенный вектор "вершина")
//    При этом (A - λI) * v_1 = 0, то есть (A - λI)^p * v_p = 0, но (A - λI)^(p-1) * v_p ≠ 0.
//    Такая цепочка соответствует жордановой клетке размера p x p с λ на диагонали
//    и единицами над ней. Столбцы матрицы S, соответствующие этой клетке, будут
//    векторы v_1, v_2, ..., v_p (именно в этом порядке).
//
//    Визуализация цепочки v_1, v_2, ..., v_p и оператора B = A - λI:
//
//        0  <--B-- v_1  <--B-- v_2  <--B-- ... <--B-- v_{p-1} <--B-- v_p
//
//    v_1 является собственным вектором: B * v_1 = 0.
//    v_p является "вершиной": B^p * v_p = 0, но B^(p-1) * v_p = v_1 ≠ 0.
//
// 3. ПОСТРОЕНИЕ БАЗИСА ДЛЯ K_λ:
//    Обозначим оператор B = A - λI.
//    Определяем индекс p_max – максимальную длину жордановой цепочки для λ.
//    Это наименьшее p, для которого dim(Ker(B^p)) достигает алгебраической кратности λ,
//    или, что то же самое, наименьшее p, для которого Ker(B^p) = Ker(B^(p+1)).
//
//    Итеративно строим цепочки, начиная с самых длинных:
//    - Для k от p_max до 1 (длина текущей цепочки):
//      - Ищем векторы v ("вершины цепочек") такие, что:
//        a) v ∈ Ker(B^k) (то есть B^k * v = 0)
//        b) v ∉ Ker(B^(k-1)) (то есть B^(k-1) * v ≠ 0). Такие векторы могут породить цепочку длины k.
//           Эти векторы выбираются из базиса Ker(B^k) / Ker(B^(k-1)) (дополняем базис Ker(B^(k-1)) до базиса Ker(B^k)).
//        c) Построенная из v цепочка {B^(k-1)v, ..., Bv, v} линейно независима со ВСЕМИ
//           уже найденными векторами для матрицы S (включая векторы из цепочек для
//           других λ и другие цепочки для текущего λ).
//      - Если такой вектор v найден:
//        - Добавляем векторы цепочки B^(k-1)v, B^(k-2)v, ..., Bv, v в матрицу S.
//        - Добавляем соответствующую жорданову клетку k x k в матрицу J.
//        - Уменьшаем количество "требуемых" векторов для данного λ на k.
//      - Повторяем, пока не наберем k_λ векторов для данного λ.
//
// ТЕХНИЧЕСКИЕ ДЕТАЛИ:
// - Решение систем линейных уравнений (для поиска Ker(B^k)): Используется метод Гаусса
//   для приведения матрицы к ступенчатому виду (RREF), из которого затем
//   извлекается базис нуль-пространства.
// - Проверка линейной независимости: Набор векторов линейно независим, если ранг
//   матрицы, составленной из этих векторов-столбцов, равен количеству векторов.
//   Ранг также находится методом Гаусса.
// - Работа с числами с плавающей точкой: Используется `long double` для точности
//   и `EPSILON` для сравнения чисел (например, для проверки, равен ли вектор нулю).
//
// СТРУКТУРА КОДА:
// - Вспомогательные функции для матриц и векторов (печать, операции).
// - Метод Гаусса и функции для нахождения ранга и базиса нуль-пространства.
// - Функция для проверки линейной независимости.
// - Функция `get_operator_index` для определения максимальной длины цепочки.
// - Основная функция `jordan_normal_form`.
// - Функция `main` для ввода данных и вывода результатов.
//
// =================================================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip> // Для std::setprecision
#include <utility> // Для std::pair
#include <string>  // Для std::string в print_matrix/vector

// Определяем константу EPSILON для сравнения чисел с плавающей точкой.
// Это маленькое число, которое используется как порог: если абсолютное значение
// разности двух чисел меньше EPSILON, они считаются равными (L для long double).
const long double EPSILON = 1e-9L;

// Псевдонимы типов для удобства: Matrix - это вектор векторов, Vector - вектор.
using Matrix = std::vector<std::vector<long double>>;
using Vector = std::vector<long double>;



// ============ Вспомогательные функции для работы с матрицами и векторами ============

/**
 * @brief Печатает матрицу в консоль.
 * @param name Имя матрицы для вывода.
 * @param m Матрица для печати.
 */
void print_matrix(const std::string& name, const Matrix& m) {
    std::cout << "Matrix " << name << ":" << std::endl;
    if (m.empty()) {
        std::cout << "  (empty)" << std::endl;
        return;
    }
    for (const auto& row : m) {
        std::cout << "  [";
        for (size_t j = 0; j < row.size(); ++j) {
            std::cout << row[j] << (j == row.size() - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
}

/**
 * @brief Печатает вектор в консоль.
 * @param name Имя вектора для вывода.
 * @param v Вектор для печати.
 */
void print_vector(const std::string& name, const Vector& v) {
    std::cout << "Vector " << name << ":" << std::endl;
    if (v.empty()) {
        std::cout << "  (empty)" << std::endl;
        return;
    }
    std::cout << "  [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << v[i] << (i == v.size() - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl << std::endl;
}

/**
 * @brief Проверяет, является ли вектор нулевым (все компоненты близки к нулю).
 * @param v Вектор для проверки.
 * @param tol Допуск (по умолчанию EPSILON).
 * @return true, если вектор нулевой, иначе false.
 */
bool is_zero_vector(const Vector& v, long double tol = EPSILON) {
    for (long double val : v) {
        if (std::abs(val) > tol) { // Если хотя бы одна компонента существенно не равна нулю
            return false;
        }
    }
    return true; // Все компоненты близки к нулю
}

/**
 * @brief Создает единичную матрицу заданного размера.
 *        Единичная матрица E (или I) - квадратная матрица, у которой на главной
 *        диагонали стоят единицы, а все остальные элементы - нули.
 * @param n Размер матрицы (n x n).
 * @return Единичная матрица размера n x n.
 */
Matrix identity_matrix(int n) {
    Matrix I(n, Vector(n, 0.0L)); // Создаем матрицу n x n, заполненную нулями
    for (int i = 0; i < n; ++i) {
        I[i][i] = 1.0L; // Устанавливаем 1 на главной диагонали
    }
    return I;
}

/**
 * @brief Умножает матрицу на скаляр.
 *        Каждый элемент матрицы умножается на заданное число.
 * @param m Исходная матрица.
 * @param scalar Скаляр для умножения.
 * @return Новая матрица, равная m * scalar.
 */
Matrix scalar_multiply_matrix(const Matrix& m, long double scalar) {
    if (m.empty()) return {};
    Matrix result = m; // Копируем исходную матрицу
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[0].size(); ++j) {
            result[i][j] *= scalar; // Умножаем каждый элемент
        }
    }
    return result;
}

/**
 * @brief Складывает две матрицы.
 *        Матрицы должны иметь одинаковые размеры.
 *        Результат C = A + B, где C_ij = A_ij + B_ij.
 * @param a Первая матрица.
 * @param b Вторая матрица.
 * @return Сумма матриц a и b.
 * @throw std::invalid_argument если размеры матриц не совпадают.
 */
Matrix matrix_add(const Matrix& a, const Matrix& b) {
    if (a.empty() || b.empty() || a.size() != b.size() || a[0].size() != b[0].size()) {
        throw std::invalid_argument("Matrix dimensions must match for addition.");
    }
    Matrix result = a; // Копируем первую матрицу
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a[0].size(); ++j) {
            result[i][j] += b[i][j]; // Поэлементно складываем
        }
    }
    return result;
}

/**
 * @brief Вычитает одну матрицу из другой.
 *        Реализуется как A - B = A + (-1 * B).
 * @param a Матрица, из которой вычитают.
 * @param b Матрица, которую вычитают.
 * @return Разность матриц a и b.
 */
Matrix matrix_subtract(const Matrix& a, const Matrix& b) {
    return matrix_add(a, scalar_multiply_matrix(b, -1.0L));
}

/**
 * @brief Умножает две матрицы.
 *        Пусть A - матрица m x k, B - матрица k x n.
 *        Результат C = A * B будет матрицей m x n.
 *        C_ij = Σ (A_il * B_lj) для l от 0 до k-1.
 * @param a Первая матрица (левый множитель).
 * @param b Вторая матрица (правый множитель).
 * @return Произведение матриц a и b.
 * @throw std::invalid_argument если размеры матриц несовместимы для умножения.
 */
Matrix matrix_multiply(const Matrix& a, const Matrix& b) {
    if (a.empty() || b.empty()) { // Проверка на пустые матрицы в начале
         if (a.empty() && b.empty()) return {};
         if (a.empty()) throw std::invalid_argument("First matrix in multiplication is empty.");
         if (b.empty()) throw std::invalid_argument("Second matrix in multiplication is empty.");
    }
    if (a[0].size() != b.size()) { // Основная проверка совместимости
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication.");
    }
    size_t m = a.size();            // Количество строк в A (и в результате)
    size_t k_common = a[0].size();  // Общая размерность: столбцы A / строки B
    size_t n_cols_b = b[0].size();  // Количество столбцов в B (и в результате)
    
    Matrix result(m, Vector(n_cols_b, 0.0L)); // Результирующая матрица m x n_cols_b, заполненная нулями
    
    for (size_t i = 0; i < m; ++i) {                // Итерация по строкам результата
        for (size_t j = 0; j < n_cols_b; ++j) {     // Итерация по столбцам результата
            for (size_t l = 0; l < k_common; ++l) { // Вычисление скалярного произведения i-й строки A на j-й столбец B
                result[i][j] += a[i][l] * b[l][j];
            }
        }
    }
    return result;
}

/**
 * @brief Умножает матрицу на вектор-столбец.
 *        Пусть M - матрица m x n, v - вектор-столбец n x 1.
 *        Результат w = M * v будет вектором-столбцом m x 1.
 *        w_i = Σ (M_ij * v_j) для j от 0 до n-1.
 * @param m Матрица.
 * @param v Вектор.
 * @return Вектор-результат умножения.
 * @throw std::invalid_argument если размеры несовместимы.
 */
Vector matrix_vector_multiply(const Matrix& m, const Vector& v) {
    if (m.empty()) {
        if (v.empty()) return {}; // Произведение пустой матрицы на пустой вектор = пустой вектор
        throw std::invalid_argument("Matrix is empty for matrix-vector multiplication.");
    }
    if (m[0].size() != v.size()) {
        throw std::invalid_argument("Matrix and vector dimensions incompatible for multiplication.");
    }
    size_t rows = m.size();
    size_t cols = m[0].size(); // == v.size()
    Vector result(rows, 0.0L); // Результирующий вектор, заполненный нулями
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[i] += m[i][j] * v[j]; // Скалярное произведение i-й строки матрицы m на вектор v
        }
    }
    return result;
}

/**
 * @brief Возводит квадратную матрицу в неотрицательную целую степень.
 *        Используется алгоритм бинарного возведения в степень (exponentiation by squaring)
 *        для эффективности.
 *        A^0 = I (единичная матрица)
 *        A^n = A * A^(n-1)
 * @param base Исходная квадратная матрица.
 * @param exp Степень (неотрицательное целое число).
 * @return Матрица base, возведенная в степень exp.
 * @throw std::invalid_argument если матрица не квадратная или степень отрицательная.
 */
Matrix matrix_power(Matrix base, int exp) { // base копируется, т.к. будет изменяться внутри (или current_power)
    if (base.empty() || base.size() != base[0].size()) {
        throw std::invalid_argument("Matrix must be square for exponentiation.");
    }
    if (exp < 0) {
        throw std::invalid_argument("Negative exponent not supported for matrix power.");
    }
    int n_dim = base.size();
    Matrix result = identity_matrix(n_dim); // A^0 = I
    if (exp == 0) return result;

    // Алгоритм бинарного возведения в степень
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = matrix_multiply(result, base);
        }
        base = matrix_multiply(base, base);
        exp /= 2;
    }
    return result;
}



// ============ Метод Гаусса и связанные с ним функции для анализа систем линейных уравнений ============

/**
 * @brief Приводит матрицу к улучшенному ступенчатому виду (RREF - Reduced Row Echelon Form)
 *        методом Гаусса-Жордана. Модифицирует матрицу на месте.
 * @param A Матрица, которую нужно привести к RREF. Модифицируется на месте.
 * @return Ранг матрицы (количество ненулевых строк в RREF, или количество пивотов).
 */
int gaussian_elimination_rref(Matrix& A) {
    if (A.empty()) return 0;
    int rows = A.size();
    int cols = A[0].size();
    int lead = 0;  // Индекс текущего ведущего столбца (пивотного столбца)
    int rank = 0;  // Ранг матрицы, будет равен количеству найденных пивотов

    // Проходим по строкам, пока не закончатся строки или столбцы
    for (int r = 0; r < rows && lead < cols;) {

        int i = r;  // Ищем строку с ненулевым элементом в столбце 'lead', начиная со строки 'r'
        while (i < rows && std::abs(A[i][lead]) < EPSILON) {
            i++;
        }

        if (i == rows) { // Если все элементы в столбце 'lead' ниже или на A[r][lead] равны нулю
            lead++;      // Переходим к следующему столбцу (`r` не инкрементируем, т.к. для строки `r` пивот ещё не найден)
            continue;
        }

        // Переставляем строку с найденным пивотом (строка 'i') на место текущей строки 'r'
        std::swap(A[i], A[r]); 
        
        // Нормализуем пивотную строку A[r]: делим её на значение пивота A[r][lead],
        // чтобы пивотный элемент стал равен 1.
        long double pivot_val = A[r][lead];
        if (std::abs(pivot_val) < EPSILON) { 
            // Этого не должно произойти, если предыдущая проверка на `i == rows` верна,
            // но на всякий случай для избежания деления на ноль.
            lead++; 
            continue;
        }
        for (int j = lead; j < cols; ++j) { // Делим все элементы строки, начиная с пивота
            A[r][j] /= pivot_val;
        }
        A[r][lead] = 1.0L; // Устанавливаем пивот точно в 1.0 для избежания ошибок округления

        // Обнуляем все остальные элементы в столбце 'lead' (кроме самого пивота)
        for (int k = 0; k < rows; ++k) {
            if (k != r) { // Для всех строк, кроме текущей пивотной строки 'r'
                long double factor = A[k][lead];     // Коэффициент, на который нужно умножить пивотную строку
                for (int j = lead; j < cols; ++j) {  // Вычитаем (factor * пивотная_строка) из строки 'k'
                    A[k][j] -= factor * A[r][j];
                }
                 A[k][lead] = 0.0L;
            }
        }
        lead++; // Переходим к следующему ведущему столбцу
        r++;    // Переходим к следующей строке для поиска пивота
        rank++; // Увеличиваем ранг, так как нашли ещё один пивот
    }
    return rank;
}

/**
 * @brief Находит базис нуль-пространства (ядра) матрицы по её RREF-форме.
 *        Нуль-пространство (ядро) Ker(M) для матрицы M - это множество всех векторов x,
 *        таких что M*x = 0.
 *
 *        Как найти базис из RREF:
 *        1. Определить пивотные (базисные) и свободные переменные.
 *           - Столбцы с пивотами соответствуют базисным переменным.
 *           - Столбцы без пивотов соответствуют свободным переменным.
 *        2. Для каждой свободной переменной (столбец `j` без пивота) построить один базисный вектор нуль-пространства `bv`:
 *           - Установить компоненту `bv[j]` (соответствующую этой свободной переменной) равной 1.
 *           - Все остальные компоненты `bv`, соответствующие *другим* свободным переменным, равны 0.
 *           - Для каждой пивотной переменной (столбец `p` с пивотом в строке `r`):
 *             `bv[p] = -rref_A[r][j]`. Это следует из уравнения `rref_A[r][p]*x_p + ... + rref_A[r][j]*x_j + ... = 0`,
 *             где `rref_A[r][p]=1`, `x_j=1`, а остальные свободные переменные 0.
 *             Тогда `x_p + rref_A[r][j]*1 = 0 => x_p = -rref_A[r][j]`.
 *
 * @param rref_A Матрица в форме RREF.
 * @return Вектор векторов, где каждый внутренний вектор - это базисный вектор нуль-пространства.
 */
std::vector<Vector> get_nullspace_basis_from_rref(const Matrix& rref_A) {
    if (rref_A.empty()) return {};
    int rows = rref_A.size();
    int cols = rref_A[0].size();
    std::vector<Vector> basis; // Сюда будем складывать базисные векторы нуль-пространства

    // Определяем, какие столбцы являются пивотными
    // pivot_col_for_row[r] = c  означает, что в строке r пивот находится в столбце c
    // is_pivot_column[c] = true означает, что столбец c является пивотным
    std::vector<int> pivot_col_for_row(rows, -1); 
    std::vector<bool> is_pivot_column(cols, false); 
    
    int current_row = 0;
    for (int c = 0; c < cols && current_row < rows; ++c) {
        if (std::abs(rref_A[current_row][c] - 1.0L) < EPSILON) { // Нашли пивот (единицу)
            // Дополнительно убедимся, что это действительно пивот (остальные в столбце нули)
            // Это уже должно быть так после gaussian_elimination_rref
            bool is_true_pivot = true;
            for(int check_r = 0; check_r < rows; ++check_r) {
                if (check_r != current_row && std::abs(rref_A[check_r][c]) > EPSILON) {
                    is_true_pivot = false;
                    break;
                }
            }
            if (is_true_pivot) {
                pivot_col_for_row[current_row] = c;
                is_pivot_column[c] = true;
                current_row++; // Переходим к следующей строке для поиска следующего пивота
            }
        }
        // Если rref_A[current_row][c] не 1, но current_row < rows, то этот столбец не пивотный
        // для текущей строки current_row. Пивот для current_row (если есть) будет правее.
        // Если rref_A[current_row][c] просто 0, то тоже не пивот.
    }
    
    // Для каждого свободного столбца (не пивотного) строим базисный вектор
    for (int j = 0; j < cols; ++j) { // j - индекс текущего столбца
        if (!is_pivot_column[j]) { // Если столбец j - свободный
            Vector basis_vector(cols, 0.0L); // Создаем новый базисный вектор
            basis_vector[j] = 1.0L; // Устанавливаем j-ю (свободную) компоненту в 1

            // Вычисляем значения для пивотных компонент этого базисного вектора
            for (int r = 0; r < rows; ++r) { // r - номер строки в RREF
                if (pivot_col_for_row[r] != -1) { // Если в строке r есть пивот
                    int pivot_c = pivot_col_for_row[r]; // Индекс пивотного столбца для строки r
                    // Из уравнения RREF[r][pivot_c]*x_pivot_c + ... + RREF[r][j]*x_j + ... = 0
                    // где x_pivot_c - пивотная переменная, x_j - свободная (=1)
                    // RREF[r][pivot_c] = 1. Остальные свободные переменные = 0.
                    // Получаем: 1*x_pivot_c + RREF[r][j]*1 = 0  => x_pivot_c = -RREF[r][j]
                    basis_vector[pivot_c] = -rref_A[r][j];
                }
            }
            basis.push_back(basis_vector); // Добавляем построенный вектор в базис
        }
    }
    return basis;
}

/**
 * @brief Находит базис нуль-пространства (ядра) для произвольной матрицы A.
 *        Сначала приводит A к RREF, затем вызывает get_nullspace_basis_from_rref.
 * @param A Исходная матрица (копируется, так как RREF её изменяет).
 * @return Базис нуль-пространства матрицы A.
 */
std::vector<Vector> get_nullspace_basis(Matrix A) { // A копируется
    if (A.empty()) return {};
    gaussian_elimination_rref(A); // A теперь в RREF-форме
    return get_nullspace_basis_from_rref(A);
}

/**
 * @brief Проверяет, является ли заданный набор векторов линейно независимым.
 *        Комбинирует уже существующие столбцы `existing_cols` с новыми `new_cols`.
 *        Метод:
 *        1. Составить матрицу M, где столбцы - это все переданные векторы.
 *        2. Если ранг матрицы M равен количеству векторов, то векторы линейно независимы.
 *           В противном случае - линейно зависимы.
 * @param existing_cols Векторы, которые уже считаются частью базиса/набора.
 * @param new_cols Новые векторы, которые нужно проверить на ЛНЗ с existing_cols.
 * @return true, если ОБЪЕДИНЕННЫЙ набор векторов линейно независим, иначе false.
 */
bool are_vectors_linearly_independent(const std::vector<Vector>& existing_cols, const std::vector<Vector>& new_cols) {
    std::vector<Vector> all_cols = existing_cols;                                           // Копируем существующие
    all_cols.insert(all_cols.end(), new_cols.begin(), new_cols.end()); // Добавляем новые

    if (all_cols.empty()) return true; // Пустой набор векторов ЛНЗ
    if (all_cols[0].empty()) { // Если первый вектор пуст (0 компонент) --> если все векторы пустые, можно считать ЛНЗ (тривиально)
        for(const auto& v : all_cols) if (!v.empty()) throw std::runtime_error("LI check: mix of empty/non-empty vectors");
        return true; 
    }

    int num_vectors = all_cols.size();    // Количество векторов в объединенном наборе
    int dim = all_cols[0].size();         // Размерность
    if (num_vectors > dim) return false;  // Если векторов больше, чем их размерность, они точно линейно зависимы
    if (num_vectors == 0) return true;    // Если после объединения нет векторов - ЛНЗ

    // Составляем матрицу M, где столбцы - это векторы из all_cols.
    // Матрица M будет иметь `dim` строк и `num_vectors` столбцов.
    Matrix M(dim, Vector(num_vectors));
    for (int j = 0; j < num_vectors; ++j) {      // j - индекс столбца (вектора)
        if (all_cols[j].size() != (size_t)dim) { // Доп. проверка на консистентность размерности
            throw std::runtime_error("Inconsistent vector dimensions for LI check.");
        }
        for (int i = 0; i < dim; ++i) { // i - индекс строки (компоненты вектора)
            M[i][j] = all_cols[j][i];
        }
    }
    
    // Приводим M к RREF и получаем ранг. M будет изменена.
    int rank = gaussian_elimination_rref(M); 
    // Векторы линейно независимы, если ранг матрицы равен количеству векторов.
    return rank == num_vectors;
}



// ============ Специальные вспомогательные функции для Жордановой формы ============

/**
 * @brief Определяет индекс оператора B = (A - λI) для собственного значения λ.
 *        Индекс p - это наименьшее натуральное число, такое что
 *        Ker(B^p) = Ker(B^(p+1)).
 *        Это также означает, что dim(Ker(B^p)) равна алгебраической кратности λ,
 *        если Ker(B^p) является корневым подпространством для λ.
 *        Индекс p соответствует максимальной длине жордановой цепочки для λ.
 *
 *        Алгоритм:
 *        Последовательно вычисляем B, B^2, B^3, ...
 *        Для каждой степени B^i находим размерность ядра dim(Ker(B^i)).
 *        Процесс останавливается, когда dim(Ker(B^i)) перестает увеличиваться
 *        или когда dim(Ker(B^i)) достигает алгебраической кратности λ.
 *
 * @param B_op Матрица оператора (A - λI).
 * @param algebraic_multiplicity Алгебраическая кратность собственного значения λ.
 * @param n_dim Размерность пространства (размер матрицы A).
 * @return Индекс p_max (максимальная длина жордановой цепочки для λ).
 */
int get_operator_index(const Matrix& B_op, int algebraic_multiplicity, int n_dim) {
    if (algebraic_multiplicity == 0) return 0; // Если кратность 0, то и индекса нет (или он 0)
    
    // Ker(B^0) = Ker(I) = {0}, поэтому dim(Ker(B^0)) = 0.
    int prev_nullity = 0; // Nullity(B^(p-1)) = n_dim - Rank(B^(p-1))

    // Итерируем по степеням p от 1. Теоретически, p не превысит algebraic_multiplicity (или n_dim).
    for (int p = 1; p <= n_dim + 1 ; ++p) { // p - текущая степень B_op
        Matrix B_p = matrix_power(B_op, p); // Вычисляем B_op^p
        Matrix rref_B_p = B_p; // Копируем для Гаусса
        int rank_B_p = gaussian_elimination_rref(rref_B_p); // Находим ранг B_op^p
        int current_nullity = n_dim - rank_B_p; // Nullity(B_op^p) = dim(Ker(B_op^p))

        // Условие стабилизации: ядро перестало расширяться
        if (current_nullity == prev_nullity) {
            // Это означает, что Ker(B^p) = Ker(B^(p-1)).
            // Следовательно, индекс равен p-1.
            return p - 1;
        }
        // Условие достижения алгебраической кратности (если она известна точно для корневого подпространства)
        // Это условие аналогично предыдущему... но на всякий случай и его проверяем:
        if (current_nullity == algebraic_multiplicity) {
            return p;
        }
        
        prev_nullity = current_nullity; // Обновляем предыдущую размерность ядра

        if (p > n_dim) { // Защита от бесконечного цикла (теоретически p <= n_dim)
            std::cerr << "Warning: Operator index search for B_op (lambda related) exceeded n=" << n_dim
                      << ". Using p=" << (p-1) << " based on last stable or n_dim." << std::endl;
            return (p-1 > 0 && p-1 <= n_dim) ? (p-1) : n_dim; // Возвращаем предыдущий p или n_dim
        }
    }
    // Если цикл завершился (например, p дошло до algebraic_multiplicity + 1 или n_dim + 1),
    // и условия выхода не сработали, это может указывать на проблему
    // (например, неверная algebraic_multiplicity).
    // Возвращаем algebraic_multiplicity как оценку "сверху", или последнюю стабильную prev_nullity.
    std::cerr << "Warning: Operator index search did not conclusively find p. "
              << "Using algebraic_multiplicity = " << algebraic_multiplicity
              << " or last known dimension of kernel for p." << std::endl;
    return algebraic_multiplicity; // Фолбэк, или можно p-1 из цикла выше
}



// ============ Основная функция для построения Жордановой нормальной формы ============

/**
 * @brief Находит Жорданову нормальную форму J и матрицу перехода S для матрицы A.
 *
 *        Алгоритм описан в начале файла. Ключевые шаги для каждого собственного значения λ:
 *        1. Формируем оператор B = A - λI.
 *        2. Находим индекс p_max = get_operator_index(B, k_λ, n) - макс. длина цепочки.
 *        3. Итерация по длинам цепочек k от p_max до 1:
 *           a. Найти базис пространства Ker(B^k).
 *           b. Для каждого вектора `v_top` из этого базиса:
 *              i.  Проверить, что `v_top` НЕ принадлежит Ker(B^(k-1)) (т.е. `B^(k-1) * v_top != 0`).
 *                  Это гарантирует, что `v_top` может быть "вершиной" цепочки длины `k`.
 *              ii. Сформировать кандидатов в цепочку: `chain = {B^(k-1)v_top, ..., Bv_top, v_top}`.
 *              iii.Проверить, что `chain` линейно независима со всеми *уже добавленными*
 *                  векторами в `S_columns_list`.
 *              iv. Если оба условия выполнены:
 *                  - Добавить векторы `chain` в `S_columns_list` (это столбцы будущей S).
 *                  - Сформировать соответствующий жорданов блок в `J_matrix`.
 *                  - Обновить счетчик найденных векторов для данного λ.
 *                  - Перейти к следующей цепочке (возможно, той же длины, если их несколько).
 *           c. Повторять, пока не будет найдено `k_λ` векторов для собственного значения λ.
 *
 * @param A Исходная квадратная матрица.
 * @param eigenvalues_list Вектор пар {собственное_значение, его_алгебраическая_кратность}.
 *                         Предполагается, что сумма кратностей равна n.
 * @return Пара {J_matrix, S_matrix}, где J - ЖНФ, S - матрица перехода.
 */
std::pair<Matrix, Matrix> jordan_normal_form(
    const Matrix& A,
    const std::vector<std::pair<long double, int>>& eigenvalues_list) { 
    
    if (A.empty()) return {{}, {}};
    int n = A.size();
    if (n == 0 || A[0].size() != (size_t)n) {
        throw std::invalid_argument("Input matrix A must be square.");
    }

    // Инициализация матриц J и S нулями
    Matrix S_matrix(n, Vector(n, 0.0L));
    Matrix J_matrix(n, Vector(n, 0.0L));
    
    // В S_columns_list будем хранить найденные векторы жорданова базиса (будущие столбцы S)
    // Порядок важен: векторы одной цепочки должны идти подряд.
    std::vector<Vector> S_columns_list;
    S_columns_list.reserve(n); // Резервируем место для n векторов

    // Индекс, с которого начинается следующий диагональный блок в J_matrix
    int J_current_block_start_idx = 0;

    // Сортируем собственные значения для консистентного вывода J и S.
    // Не обязательно для корректности, но полезно для воспроизводимости.
    std::vector<std::pair<long double, int>> sorted_eigenvalues = eigenvalues_list;
    std::sort(sorted_eigenvalues.begin(), sorted_eigenvalues.end(), 
        [](const auto&a, const auto&b){
            // Сортируем по значению собственного числа
            if (std::abs(a.first - b.first) > EPSILON * 100) { // Увеличенный EPSILON для сравнения ключей сортировки
                 return a.first < b.first;
            }
            // Для одинаковых собственных чисел (если бы это было возможно при правильной группировке)
            // можно было бы сортировать по кратности, но это не должно влиять.
            return a.second > b.second; 
        }
    );

    // --- Основной цикл по каждому УНИКАЛЬНОМУ собственному значению ---
    for (const auto& eig_pair : sorted_eigenvalues) {
        long double lambda = eig_pair.first;  // Текущее собственное значение
        int multiplicity = eig_pair.second;   // Его алгебраическая кратность

        if (multiplicity == 0) continue; // Пропускаем, если кратность 0

        // Формируем оператор B = A - λI
        Matrix B_operator = matrix_subtract(A, scalar_multiply_matrix(identity_matrix(n), lambda));
        
        // Находим индекс p_max - максимальную длину жордановой цепочки для этого λ
        int p_max_chain_length = get_operator_index(B_operator, multiplicity, n);
        
        // Если p_max_chain_length = 0, а multiplicity > 0, это странная ситуация.
        // dim(Ker(B^0)) = 0, не может быть равно multiplicity. Значит, get_operator_index вернул 0.
        // Это может случиться, если algebraic_multiplicity в get_operator_index было 0.
        // Но мы отфильтровали multiplicity=0.
        // Если B_operator -- нулевая матрица, то Ker(B_operator) = R^n.
        // Тогда p_max_chain_length будет 1.
        if (p_max_chain_length == 0 && multiplicity > 0) {
            std::cerr << "Warning or Info: For lambda=" << lambda << " (multiplicity " << multiplicity
                      << "), max chain length (p_index) is 0. "
                      << "This implies B=(A-lambda*I) might be special (e.g., zero on the relevant subspace "
                      << "and A is already diagonal for this lambda, or multiplicity is actually 0)." << std::endl;
            // Если A было lambda*I, то B=0, Ker(B)=R^n, p_max=1. Get_operator_index должен это корректно обработать.
            // Если p_max_chain_length тут 0, возможно, get_operator_index требует доработки для такого случая
            // или исходные данные некорректны (multiplicity > 0, но лямбда не собственное число).
            // Для простоты, если p_max=0, а кратность >0, то ставим p_max=1, т.к. хотя бы один собств. вектор есть.
             if (multiplicity > 0) p_max_chain_length = 1; 
        }


        int vecs_found_for_this_lambda = 0; // Счетчик найденных базисных векторов для текущего λ

        // --- Цикл по возможным длинам цепочек k (от максимальной p_max_chain_length до 1) ---
        for (int k_chain_length = p_max_chain_length; k_chain_length >= 1; --k_chain_length) { 
            if (vecs_found_for_this_lambda == multiplicity) break; // Нашли достаточно векторов для этого λ

            // Нам нужны векторы v из Ker(B^k) такие, что v НЕ из Ker(B^(k-1)).
            // То есть, B^k * v = 0, но B^(k-1) * v != 0.
            
            Matrix B_pow_k = matrix_power(B_operator, k_chain_length);
            Matrix B_pow_k_minus_1 = (k_chain_length == 0) ? identity_matrix(n) // B^0 = I
                                         : matrix_power(B_operator, k_chain_length - 1);

            // Находим полный базис для Ker(B^k)
            std::vector<Vector> basis_ker_Bk = get_nullspace_basis(B_pow_k);
             if (basis_ker_Bk.empty() && k_chain_length > 0 && multiplicity > vecs_found_for_this_lambda) {
                  // Если ядро B_pow_k нулевое, а мы ожидаем векторы, что-то не так.
                  // Либо p_max_chain_length был вычислен неверно, либо ошибка в get_nullspace_basis.
                  std::cerr << "Warning: Ker(B^" << k_chain_length << ") for lambda " << lambda 
                            << " is empty or not found, but more vectors expected." << std::endl;
                  continue; 
             }

            // --- Итерация по кандидатам в "вершины" цепочек из базиса Ker(B^k) ---
            for (const auto& v_candidate_top : basis_ker_Bk) {
                if (is_zero_vector(v_candidate_top, EPSILON * 10)) continue; // Пропускаем нулевой вектор, если он попал в базис
                if (vecs_found_for_this_lambda == multiplicity) break; // Уже нашли все для этого λ

                // 1. Проверяем, что v_candidate_top НЕ лежит в Ker(B^(k-1))
                //    То есть, B^(k-1) * v_candidate_top должно быть НЕ равно нулю.
                //    Этот вектор B^(k-1) * v_candidate_top будет первым вектором в цепочке (собственным вектором).
                Vector first_vec_in_chain = matrix_vector_multiply(B_pow_k_minus_1, v_candidate_top);
                
                if (is_zero_vector(first_vec_in_chain, EPSILON)) {
                    // Этот v_candidate_top лежит в Ker(B^(k-1)), значит, он не может начать
                    // цепочку длины k_chain_length. Он может начать более короткую цепочку.
                    // (Это будет обработано на следующих итерациях по k_chain_length).
                    continue;
                }

                // 2. Формируем саму цепочку, начиная с "вершины" v_candidate_top
                //    Цепочка для S будет: { B^(k-1)v, B^(k-2)v, ..., Bv, v }
                //    B^0*v = v
                //    B^1*v
                //    ...
                //    B^(k-1)*v
                std::vector<Vector> current_chain_vectors;
                current_chain_vectors.reserve(k_chain_length);
                for (int j = 0; j < k_chain_length; ++j) { 
                    // Вектор v_i = B^( (k-1) - (i-1) ) * v_top  для i=1..k (индекс в цепочке)
                    // В S_columns_list они должны идти в порядке v_1, v_2, ... v_k
                    // j=0 -> B_pow = k-1.   v_1 = B^(k-1)v_top
                    // j=1 -> B_pow = k-2.   v_2 = B^(k-2)v_top
                    // ...
                    // j=k-1 -> B_pow = 0.   v_k = B^0 v_top = v_top
                    current_chain_vectors.push_back(matrix_vector_multiply(matrix_power(B_operator, k_chain_length - 1 - j), v_candidate_top));
                }
                
                // 3. Проверяем линейную независимость новой цепочки со ВСЕМИ УЖЕ НАЙДЕННЫМИ
                //    векторами в S_columns_list.
                if (are_vectors_linearly_independent(S_columns_list, current_chain_vectors)) {
                    // Если ЛНЗ, то цепочка подходит!
                    // Добавляем векторы этой цепочки в общий список S_columns_list
                    S_columns_list.insert(S_columns_list.end(), current_chain_vectors.begin(), current_chain_vectors.end());
                    
                    // Формируем соответствующий жорданов блок k_chain_length x k_chain_length в J_matrix
                    // Он начинается с J_matrix[J_current_block_start_idx][J_current_block_start_idx]
                    for (int i = 0; i < k_chain_length; ++i) {
                        J_matrix[J_current_block_start_idx + i][J_current_block_start_idx + i] = lambda; // λ на диагонали
                        if (i < k_chain_length - 1) {
                            // Единицы над главной диагональю внутри блока
                            J_matrix[J_current_block_start_idx + i][J_current_block_start_idx + i + 1] = 1.0L;
                        }
                    }
                    J_current_block_start_idx += k_chain_length; // Сдвигаем начало следующего блока
                    vecs_found_for_this_lambda += k_chain_length; // Обновляем счетчик найденных векторов
                }
                // Если не ЛНЗ, то этот v_candidate_top (или порожденная им цепочка) является
                // линейной комбинацией уже найденных векторов и/или других векторов из Ker(B^k).
                // Ищем другой v_candidate_top.
            }
        } // конец цикла по k_chain_length

        // Проверка, все ли векторы для данного λ найдены
        if (vecs_found_for_this_lambda != multiplicity) {
             std::cerr << "Warning: For eigenvalue " << lambda 
                       << " (multiplicity " << multiplicity << "), "
                       << "expected " << multiplicity << " basis vectors, but found " 
                       << vecs_found_for_this_lambda << "." << std::endl;
             std::cerr << "         This might happen due to numerical precision issues or "
                       << "incorrect input eigenvalues/multiplicities." << std::endl;
        }
    } // конец цикла по собственным значениям (eig_pair)

    // Финальная проверка общего количества найденных векторов
    if (S_columns_list.size() != (size_t)n) {
        std::cerr << "Error: Expected to find " << n << " basis vectors for S, but found " 
                  << S_columns_list.size() << "." << std::endl;
        std::cerr << "         The S matrix will be incomplete or incorrect. "
                  << "The J matrix might also be incorrect." << std::endl;
        // Можно либо бросить исключение, либо вернуть "как есть" с неполной S.
        // Если S_columns_list.size() < n, то S_matrix будет иметь нулевые столбцы в конце,
        // что сделает S сингулярной.
    } else {
          std::cout << "Successfully found " << n << " basis vectors for S matrix." << std::endl;
    }

    // Заполняем матрицу S_matrix из списка векторов S_columns_list
    // Каждый вектор из S_columns_list становится столбцом в S_matrix
    for (size_t j = 0; j < S_columns_list.size(); ++j) { // j - номер столбца в S_matrix
        if (S_columns_list[j].size() != (size_t)n) { // Проверка размерности вектора
             throw std::runtime_error("Internal error: S column vector has incorrect dimension during S_matrix assembly.");
        }
        for (int i = 0; i < n; ++i) { // i - номер строки в S_matrix
            if (j < (size_t)n) { // Доп. предохранитель, чтобы не выйти за границы S_matrix, если S_columns_list < n
                 S_matrix[i][j] = S_columns_list[j][i];
            }
        }
    }
    
    return {J_matrix, S_matrix};
}



// ============ Главная функция программы ============
int main() {
    // Ускорение стандартного ввода-вывода (может быть полезно для больших матриц)
    std::ios_base::sync_with_stdio(false); 
    std::cin.tie(NULL);

    // Устанавливаем точность вывода для чисел с плавающей точкой
    std::cout << std::fixed << std::setprecision(3); 

    int n; // Размерность квадратной матрицы
    std::cout << "Enter the size of the square matrix (n):" << std::endl; 
    std::cin >> n;

    if (n <= 0) {
        std::cerr << "Error: Matrix size must be a positive integer." << std::endl;
        return 1;
    }
    if (n > 20) { // Примерное ограничение, выше может быть очень долго
        std::cout << "Warning: For n > 20, computations can be very slow due to matrix exponentiation and Gaussian elimination." << std::endl;
    }


    Matrix A_input(n, Vector(n)); // Инициализация входной матрицы A
    std::cout << "Enter the elements of the " << n << "x" << n << " matrix A (row by row, elements in a row separated by spaces):" << std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(std::cin >> A_input[i][j])) { // Чтение и проверка корректности ввода
                 std::cerr << "Error: Invalid input for matrix element." << std::endl; return 1;
            }
        }
    }

    // Считываем собственные значения. Пользователь вводит все n штук, с повторениями для кратностей.
    std::vector<long double> raw_eigenvalues_input;
    raw_eigenvalues_input.reserve(n); // Резервируем память
    std::cout << "Enter all " << n << " eigenvalues (including multiplicities, separated by spaces):" << std::endl;
    long double val;
    for (int i = 0; i < n; ++i) {
        if (!(std::cin >> val)) { // Чтение и проверка
             std::cerr << "Error: Invalid input for eigenvalue." << std::endl; return 1;
        }
        raw_eigenvalues_input.push_back(val);
    }

    if (raw_eigenvalues_input.size() != (size_t)n) { // Дополнительная проверка количества
        std::cerr << "Error: Number of eigenvalues entered (" << raw_eigenvalues_input.size() 
                  << ") does not match matrix dimension (" << n << ")." << std::endl;
        return 1;
    }
    
    // Обрабатываем сырой список собственных значений, чтобы получить уникальные значения и их кратности
    std::vector<std::pair<long double, int>> eigenvalues_with_multiplicity;
    if (!raw_eigenvalues_input.empty()) {
        // Сортируем, чтобы сгруппировать одинаковые (численно близкие) значения
        std::sort(raw_eigenvalues_input.begin(), raw_eigenvalues_input.end());

        // Первое собственное значение начинаем считать
        eigenvalues_with_multiplicity.push_back({raw_eigenvalues_input[0], 0}); // Кратность пока 0
        for (long double eig_val_from_input : raw_eigenvalues_input) {
            // Сравниваем текущее значение из сырого списка с последним уникальным, которое мы обрабатываем
            if (std::abs(eig_val_from_input - eigenvalues_with_multiplicity.back().first) < EPSILON * 10) { // Немного увеличенный EPSILON для группировки
                eigenvalues_with_multiplicity.back().second++; // Увеличиваем кратность текущего уникального с.з.
            } else {
                // Нашли новое уникальное собственное значение
                eigenvalues_with_multiplicity.push_back({eig_val_from_input, 1}); // Добавляем его с кратностью 1
            }
        }
    }
    
    // Проверяем, что сумма кратностей равна размерности матрицы
    int total_multiplicity_check = 0;
    for(const auto& p : eigenvalues_with_multiplicity) total_multiplicity_check += p.second;
    if(total_multiplicity_check != n) {
        std::cerr << "Critical Warning: Sum of multiplicities of processed eigenvalues (" << total_multiplicity_check 
                  << ") does not match matrix dimension (" << n << ")." << std::endl;
        std::cerr << "This usually indicates an error in the provided eigenvalues (e.g., one was missed or miscounted)." << std::endl;
        return 1;
    }


    std::cout << "\n--- Input Summary ---" << std::endl;
    print_matrix("A_input", A_input);
    std::cout << "Processed Eigenvalues with multiplicities (derived from input):" << std::endl;
    for(const auto& pair_val : eigenvalues_with_multiplicity) {
        std::cout << "  Lambda = " << pair_val.first << ", Multiplicity = " << pair_val.second << std::endl;
    }
    std::cout << std::endl;

    std::cout << "--- Computing Jordan Normal Form ---" << std::endl;
    try {
        // Вызываем основную функцию для получения J и S
        auto [J_computed, S_computed] = jordan_normal_form(A_input, eigenvalues_with_multiplicity);
        std::cout << "\n--- Results ---" << std::endl;
        print_matrix("J_computed (Jordan Form)", J_computed);
        print_matrix("S_computed (Transition Matrix)", S_computed);

    } catch (const std::exception& e) {
        // Ловим возможные исключения (например, из матричных операций)
        std::cerr << "An error occurred during computation: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
