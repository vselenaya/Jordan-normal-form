#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <utility>
#include <string>
#include <type_traits> 
#include <optional>

static constexpr int64_t PrimeForZpArithmetic = 239;  // Можно поменять на произвольное простое число
static constexpr int DoublePrintPrecision = 8;        // Количество знаков после запятой
static constexpr long double GLOBAL_EPSILON = 1e-8;   // Точность сравнения чисел с плавающей точкой



// Класс для модульной арифметики (можем использовать или его, или обычный вещественные числа... смысл первого в том, что все вычисления
// по простому модулю (в том числе деление) -- без потери точности!)
template <typename IntType, IntType Modulus>
class Modular {
public:
    IntType value;
    Modular(IntType val = 0) {
        value = val % Modulus;
        if (value < 0) {
            value += Modulus;
        }
    }
    explicit operator IntType() const { return value; }
    IntType get() const { return value; }
    Modular& operator+=(const Modular& other) {
        value += other.value;
        if (value >= Modulus) value -= Modulus;
        return *this;
    }
    Modular& operator-=(const Modular& other) {
        value -= other.value;
        if (value < 0) value += Modulus;
        return *this;
    }
    Modular& operator*=(const Modular& other) {
        value = (static_cast<long long>(value) * other.value) % Modulus;
        return *this;
    }
    static Modular power(Modular base, IntType exp) {
        Modular res(1);
        base.value %= Modulus;
        while (exp > 0) {
            if (exp % 2 == 1) res *= base;
            base *= base;
            exp /= 2;
        }
        return res;
    }
    Modular inverse() const {
        if (value == 0) throw std::runtime_error("Modular division by zero.");
        return power(*this, Modulus - 2);
    }
    Modular& operator/=(const Modular& other) {
        *this *= other.inverse();
        return *this;
    }
    friend Modular operator+(Modular lhs, const Modular& rhs) { return lhs += rhs; }
    friend Modular operator-(Modular lhs, const Modular& rhs) { return lhs -= rhs; }
    friend Modular operator*(Modular lhs, const Modular& rhs) { return lhs *= rhs; }
    friend Modular operator/(Modular lhs, const Modular& rhs) { return lhs /= rhs; }
    friend bool operator==(const Modular& lhs, const Modular& rhs) { return lhs.value == rhs.value; }
    friend bool operator!=(const Modular& lhs, const Modular& rhs) { return lhs.value != rhs.value; }
    friend bool operator<(const Modular& lhs, const Modular& rhs) { return lhs.value < rhs.value; }
    friend std::ostream& operator<<(std::ostream& os, const Modular& m) { os << m.value; return os; }
    friend std::istream& operator>>(std::istream& is, Modular& m) {
        IntType val_in; is >> val_in; m = Modular(val_in); return is;
    }
};

template <typename T> using Matrix = std::vector<std::vector<T>>;
template <typename T> using Vector = std::vector<T>;
template <typename T> void print_matrix(const std::string& name, const Matrix<T>& m); // Объявления
template <typename T> bool is_zero_scalar(const T& val);
template <typename T> bool are_scalars_equal(const T& a, const T& b);
template <typename T> Matrix<T> identity_matrix(int n);
template <typename T> Matrix<T> matrix_multiply(const Matrix<T>& a, const Matrix<T>& b);

// --- Функция обращения матрицы ---
/**
 * @brief Находит обратную матрицу M_inv для квадратной матрицы M.
 *        Использует метод Гаусса-Жордана на расширенной матрице [M | I].
 *        Если M вырождена, возвращает std::nullopt.
 * @tparam T Тип элементов матрицы (должен поддерживать полевые операции).
 * @param M Квадратная матрица для обращения.
 * @return std::optional<Matrix<T>>, содержащий обратную матрицу, или std::nullopt.
 */
template <typename T>
std::optional<Matrix<T>> matrix_inverse(const Matrix<T>& M) {
    if (M.empty() || M.size() != M[0].size()) {
        // Не квадратная или пустая
        // Можно бросить исключение или просто вернуть nullopt
        // throw std::invalid_argument("Matrix must be square and non-empty for inversion.");
        return std::nullopt; 
    }
    int n = M.size();

    // Создаем расширенную матрицу [M | I]
    Matrix<T> augmented_matrix(n, Vector<T>(2 * n));
    Matrix<T> I = identity_matrix<T>(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented_matrix[i][j] = M[i][j];
            augmented_matrix[i][j + n] = I[i][j];
        }
    }

    // Применяем метод Гаусса-Жордана к augmented_matrix
    // Цель: превратить левую часть (первые n столбцов) в единичную матрицу.
    // Правая часть (последние n столбцов) станет обратной матрицей.
    int lead = 0;
    for (int r = 0; r < n && lead < n; /* r, lead инкрементируются внутри */ ) {
        int i = r;
        while (i < n && is_zero_scalar(augmented_matrix[i][lead])) {
            i++;
        }

        if (i == n) { // Пивот не найден в столбце lead ниже или на строке r
            // Это означает, что левая часть не может быть приведена к единичной,
            // т.е. исходная матрица M сингулярна.
            return std::nullopt; 
        }

        std::swap(augmented_matrix[i], augmented_matrix[r]);
        
        T pivot_val = augmented_matrix[r][lead];
        // Если pivot_val == T(0) после swap, и i было < n, то это тоже проблема.
        // Но предыдущая проверка `if (i == n)` должна это покрывать для `lead < n`.
        if (is_zero_scalar(pivot_val)) return std::nullopt; // Сингулярна

        T inv_pivot_val = T(1) / pivot_val;
        for (int j = lead; j < 2 * n; ++j) { // Обрабатываем всю расширенную строку
            augmented_matrix[r][j] *= inv_pivot_val;
        }
        augmented_matrix[r][lead] = T(1); // Гарантируем единицу на пивоте

        for (int k = 0; k < n; ++k) {
            if (k != r) {
                T factor = augmented_matrix[k][lead];
                for (int j = lead; j < 2 * n; ++j) {
                    augmented_matrix[k][j] -= factor * augmented_matrix[r][j];
                }
                augmented_matrix[k][lead] = T(0); // Гарантируем нули в пивотном столбце
            }
        }
        lead++;
        r++;
    }
    
    // После цикла, если lead < n, это значит, что не все столбцы левой части стали пивотными,
    // т.е. ранг < n, матрица сингулярна. Хотя это должно было быть отловлено ранее.
     if (lead < n) { 
         // Проверим еще раз ранг левой половины.
         // Можно было бы вызвать gaussian_elimination_rref на копии левой половины M
         // и проверить ранг == n, но логика выше должна была это отловить.
         // Этот блок может быть избыточным.
         return std::nullopt;
     }

    // Извлекаем правую часть (обратную матрицу)
    Matrix<T> M_inv(n, Vector<T>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M_inv[i][j] = augmented_matrix[i][j + n];
        }
    }
    return M_inv;
}

// --- Функция верификации ---
/**
 * @brief Проверяет правильность разложения A = S * J * S_inv.
 * @tparam T Тип элементов.
 * @param A Исходная матрица.
 * @param J Жорданова форма.
 * @param S Матрица перехода.
 */
template <typename T>
void verify_jordan_decomposition(const Matrix<T>& A, const Matrix<T>& J, const Matrix<T>& S) {
    std::cout << "\n--- Verification ---" << std::endl;
    if (S.empty() || S.size() != A.size()) {
        std::cout << "S matrix is empty or has incorrect dimensions. Cannot verify." << std::endl;
        return;
    }

    std::optional<Matrix<T>> S_inv_opt = matrix_inverse<T>(S);

    if (S_inv_opt) {
        Matrix<T> S_inv = *S_inv_opt;
        std::cout << "S_inv successfully computed." << std::endl;
        // print_matrix<T>("S_inv_calc", S_inv); // Можно раскомментировать для отладки

        Matrix<T> S_J = matrix_multiply<T>(S, J);
        Matrix<T> S_J_S_inv = matrix_multiply<T>(S_J, S_inv);
        // print_matrix<T>("S_J_S_inv_calc", S_J_S_inv); // Для отладки

        bool success = true;
        if (A.size() != S_J_S_inv.size() || (!A.empty() && A[0].size() != S_J_S_inv[0].size())) {
            success = false;
            std::cout << "Dimension mismatch between A and S*J*S_inv." << std::endl;
        } else {
            for (size_t i = 0; i < A.size(); ++i) {
                for (size_t j = 0; j < A[0].size(); ++j) {
                    if (!are_scalars_equal(A[i][j], S_J_S_inv[i][j])) {
                        success = false;
                        break;
                    }
                }
                if (!success) break;
            }
        }

        if (success) {
            std::cout << "Verification SUCCESS: A approx S * J * S^-1" << std::endl;
        } else {
            std::cout << "Verification FAILED: A != S * J * S^-1" << std::endl;
            // Можно вывести разницу, если нужно
            // Matrix<T> Diff = matrix_subtract(A, S_J_S_inv);
            // print_matrix<T>("Difference A - SJS^-1", Diff);
        }

    } else {
        std::cout << "S is singular (or non-invertible). Cannot compute S_inv." << std::endl;
        std::cout << "Attempting to verify A * S == S * J ..." << std::endl;

        Matrix<T> A_S = matrix_multiply<T>(A, S);
        Matrix<T> S_J = matrix_multiply<T>(S, J);

        bool success_alt = true;
        if (A_S.size() != S_J.size() || (!A_S.empty() && A_S[0].size() != S_J[0].size())) {
            success_alt = false;
            std::cout << "Dimension mismatch between A*S and S*J." << std::endl;
        } else {
            for (size_t i = 0; i < A_S.size(); ++i) {
                for (size_t j = 0; j < A_S[0].size(); ++j) {
                    if (!are_scalars_equal(A_S[i][j], S_J[i][j])) {
                        success_alt = false;
                        break;
                    }
                }
                if (!success_alt) break;
            }
        }
        
        if (success_alt) {
            std::cout << "Verification SUCCESS (alternative form): A * S approx S * J" << std::endl;
        } else {
            std::cout << "Verification FAILED (alternative form): A * S != S * J" << std::endl;
        }
    }
}


// --- Далее копируем implementations для print_matrix и др. шаблонных функций ---
// ... Это важно, так как они используются в matrix_inverse и verify.
// --- Utility Functions (Шаблонизированные) ---
template <typename T>
void print_matrix(const std::string& name, const Matrix<T>& m) {
    std::cout << "Matrix " << name << ":" << std::endl;
    if (m.empty()) {
        std::cout << "  (empty)" << std::endl;
        return;
    }
    for (const auto& row : m) {
        std::cout << "  [";
        for (size_t j = 0; j < row.size(); ++j) {
            if constexpr (std::is_same_v<T, long double>) {
                 std::cout << std::fixed << std::setprecision(6) << row[j];
            } else {
                 std::cout << row[j];
            }
            std::cout << (j == row.size() - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
}

template <typename T>
void print_vector(const std::string& name, const Vector<T>& v) {
    std::cout << "Vector " << name << ":" << std::endl;
    if (v.empty()) {
        std::cout << "  (empty)" << std::endl;
        return;
    }
    std::cout << "  [";
    for (size_t i = 0; i < v.size(); ++i) {
         if constexpr (std::is_same_v<T, long double>) {
             std::cout << std::fixed << std::setprecision(6) << v[i];
        } else {
             std::cout << v[i];
        }
        std::cout << (i == v.size() - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl << std::endl;
}


template <typename T>
bool is_zero_scalar(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::abs(val) < GLOBAL_EPSILON;
    } else { 
        return val == T(0);
    }
}

template <typename T>
bool are_scalars_equal(const T& a, const T& b) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::abs(a - b) < GLOBAL_EPSILON;
    } else {
        return a == b;
    }
}


template <typename T>
bool is_zero_vector(const Vector<T>& v) {
    for (const T& val : v) {
        if (!is_zero_scalar(val)) {
            return false;
        }
    }
    return true;
}

template <typename T>
Matrix<T> identity_matrix(int n) {
    Matrix<T> I(n, Vector<T>(n, T(0)));
    for (int i = 0; i < n; ++i) {
        I[i][i] = T(1);
    }
    return I;
}

template <typename T>
Matrix<T> scalar_multiply_matrix(const Matrix<T>& m, T scalar) {
    if (m.empty()) return {};
    Matrix<T> result = m;
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[0].size(); ++j) {
            result[i][j] *= scalar;
        }
    }
    return result;
}

template <typename T>
Matrix<T> matrix_add(const Matrix<T>& a, const Matrix<T>& b) {
     if (a.empty() || b.empty() || a.size() != b.size() || (!a.empty() && !b.empty() && a[0].size() != b[0].size())) {
        if (a.empty() && b.empty()) return {};
        throw std::invalid_argument("Matrix dimensions must match for addition.");
    }
    Matrix<T> result = a;
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a[0].size(); ++j) {
            result[i][j] += b[i][j];
        }
    }
    return result;
}

template <typename T>
Matrix<T> matrix_subtract(const Matrix<T>& a, const Matrix<T>& b) {
    return matrix_add(a, scalar_multiply_matrix(b, T(-1)));
}

template <typename T>
Matrix<T> matrix_multiply(const Matrix<T>& a, const Matrix<T>& b) {
    if (a.empty() || b.empty()) {
         if (a.empty() && b.empty()) return {};
         if (a.empty()) throw std::invalid_argument("First matrix in multiplication is empty.");
         if (b.empty()) throw std::invalid_argument("Second matrix in multiplication is empty.");
    }
    if (a[0].size() != b.size()) {
        throw std::invalid_argument("Matrix dimensions (" + std::to_string(a.size()) + "x" + std::to_string(a[0].size()) + 
                                  " and " + std::to_string(b.size()) + "x" + std::to_string(b[0].size()) + 
                                  ") incompatible for multiplication.");
    }
    size_t m = a.size();      
    size_t k_common = a[0].size(); 
    size_t n_cols_b = b[0].size(); 
    
    Matrix<T> result(m, Vector<T>(n_cols_b, T(0)));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n_cols_b; ++j) {
            for (size_t l = 0; l < k_common; ++l) {
                result[i][j] += a[i][l] * b[l][j];
            }
        }
    }
    return result;
}

template <typename T>
Vector<T> matrix_vector_multiply(const Matrix<T>& m, const Vector<T>& v) {
    if (m.empty()) {
        if (v.empty()) return {};
        throw std::invalid_argument("Matrix is empty for matrix-vector multiplication.");
    }
    if (m[0].size() != v.size()) {
        throw std::invalid_argument("Matrix and vector dimensions incompatible for multiplication.");
    }
    size_t rows = m.size();
    size_t cols = m[0].size();
    Vector<T> result(rows, T(0));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[i] += m[i][j] * v[j];
        }
    }
    return result;
}

template <typename T>
Matrix<T> matrix_power(Matrix<T> base, int exp) {
    if (base.empty() || base.size() != base[0].size()) {
        throw std::invalid_argument("Matrix must be square for exponentiation.");
    }
    if (exp < 0) {
        throw std::invalid_argument("Negative exponent not supported for matrix power.");
    }
    int n_dim = base.size();
    Matrix<T> result = identity_matrix<T>(n_dim);
    if (exp == 0) return result;

    while (exp > 0) {
        if (exp % 2 == 1) {
            result = matrix_multiply(result, base);
        }
        base = matrix_multiply(base, base);
        exp /= 2;
    }
    return result;
}

template <typename T> int gaussian_elimination_rref(Matrix<T>& A);
template <typename T> std::vector<Vector<T>> get_nullspace_basis_from_rref(const Matrix<T>& rref_A);
template <typename T> std::vector<Vector<T>> get_nullspace_basis(Matrix<T> A);
template <typename T> bool are_vectors_linearly_independent(const std::vector<Vector<T>>& existing_cols, const std::vector<Vector<T>>& new_cols);
template <typename T> int get_operator_index(const Matrix<T>& B_op, int algebraic_multiplicity, int n_dim);
template <typename T> std::pair<Matrix<T>, Matrix<T>> jordan_normal_form(const Matrix<T>& A, const std::vector<std::pair<T, int>>& eigenvalues_list_typed);

// --- Gaussian Elimination (Шаблонизированная) ---
template <typename T>
int gaussian_elimination_rref(Matrix<T>& A) {
    if (A.empty()) return 0;
    int rows = A.size();
    int cols = A[0].size();
    int lead = 0; 
    int rank = 0;

    for (int r = 0; r < rows && lead < cols; ) {
        int i = r;
        while (i < rows && is_zero_scalar(A[i][lead])) {
            i++;
        }

        if (i == rows) { 
            lead++;      
            continue;
        }

        std::swap(A[i], A[r]); 
        
        T pivot_val = A[r][lead];
        if (is_zero_scalar(pivot_val)) {
             lead++; continue;
        }

        T inv_pivot_val = T(1) / pivot_val; 
        for (int j = lead; j < cols; ++j) { 
            A[r][j] *= inv_pivot_val; 
        }
        A[r][lead] = T(1); 
        for (int k = 0; k < rows; ++k) {
            if (k != r) { 
                T factor = A[k][lead]; 
                for (int j = lead; j < cols; ++j) {
                    A[k][j] -= factor * A[r][j];
                }
                 A[k][lead] = T(0); 
            }
        }
        lead++; 
        r++;    
        rank++; 
    }
    return rank;
}

template <typename T>
std::vector<Vector<T>> get_nullspace_basis_from_rref(const Matrix<T>& rref_A) {
    if (rref_A.empty()) return {};
    int rows = rref_A.size();
    int cols = rref_A[0].size();
    std::vector<Vector<T>> basis; 

    std::vector<int> pivot_col_for_row(rows, -1); 
    std::vector<bool> is_pivot_column(cols, false); 
    
    int current_row = 0;
    for (int c = 0; c < cols && current_row < rows; ++c) {
        if (are_scalars_equal(rref_A[current_row][c], T(1))) { 
            bool is_true_pivot = true;
            for(int check_r = 0; check_r < rows; ++check_r) {
                if (check_r != current_row && !is_zero_scalar(rref_A[check_r][c])) {
                    is_true_pivot = false; break;
                }
            }
            if (is_true_pivot) {
                pivot_col_for_row[current_row] = c;
                is_pivot_column[c] = true;
                current_row++; 
            }
        }
    }
    
    for (int j = 0; j < cols; ++j) { 
        if (!is_pivot_column[j]) { 
            Vector<T> basis_vector(cols, T(0)); 
            basis_vector[j] = T(1); 
            for (int r_idx = 0; r_idx < rows; ++r_idx) { 
                if (pivot_col_for_row[r_idx] != -1) { 
                    int pivot_c = pivot_col_for_row[r_idx]; 
                    basis_vector[pivot_c] = T(0) - rref_A[r_idx][j];
                }
            }
            basis.push_back(basis_vector); 
        }
    }
    return basis;
}

template <typename T>
std::vector<Vector<T>> get_nullspace_basis(Matrix<T> A) { 
    if (A.empty()) return {};
    gaussian_elimination_rref<T>(A); 
    return get_nullspace_basis_from_rref<T>(A);
}

template <typename T>
bool are_vectors_linearly_independent(const std::vector<Vector<T>>& existing_cols, const std::vector<Vector<T>>& new_cols) {
    std::vector<Vector<T>> all_cols = existing_cols; 
    all_cols.insert(all_cols.end(), new_cols.begin(), new_cols.end()); 

    if (all_cols.empty()) return true; 
    if (all_cols[0].empty()) { 
        for(const auto& v : all_cols) if (!v.empty()) throw std::runtime_error("LI check: mix of empty/non-empty vectors");
        return true; 
    }

    int num_vectors = all_cols.size(); 
    int dim = all_cols[0].size();      
    if (num_vectors > dim) return false;
    if (num_vectors == 0) return true; 

    Matrix<T> M(dim, Vector<T>(num_vectors));
    for (int j_vec = 0; j_vec < num_vectors; ++j_vec) { 
        if (all_cols[j_vec].size() != (size_t)dim) { 
            throw std::runtime_error("Inconsistent vector dimensions for LI check.");
        }
        for (int i_comp = 0; i_comp < dim; ++i_comp) { 
            M[i_comp][j_vec] = all_cols[j_vec][i_comp];
        }
    }
    
    int rank = gaussian_elimination_rref<T>(M); 
    return rank == num_vectors;
}

template <typename T>
int get_operator_index(const Matrix<T>& B_op, int algebraic_multiplicity, int n_dim) {
    if (algebraic_multiplicity == 0) return 0; 
    int prev_nullity = 0; 

    for (int p = 1; p <= n_dim + 1 ; ++p) { 
        Matrix<T> B_p = matrix_power<T>(B_op, p); 
        Matrix<T> rref_B_p = B_p; 
        int rank_B_p = gaussian_elimination_rref<T>(rref_B_p); 
        int current_nullity = n_dim - rank_B_p; 

        if (current_nullity == prev_nullity && p > 1) {
            return p - 1;
        }
        if (current_nullity == algebraic_multiplicity) {
            return p;
        }
        prev_nullity = current_nullity; 
        if (p > n_dim) { 
             std::cerr << "Warning: Operator index search for B_op (lambda related) exceeded n=" << n_dim
                       << ". Using p=" << (p-1) << " based on last stable or n_dim." << std::endl;
            return (p-1 > 0 && p-1 <= n_dim) ? (p-1) : n_dim; 
        }
    }
    std::cerr << "Warning: Operator index search did not conclusively find p. "
              << "Using algebraic_multiplicity = " << algebraic_multiplicity
              << " or last known dimension of kernel for p." << std::endl;
    return algebraic_multiplicity; 
}



template <typename T>
std::pair<Matrix<T>, Matrix<T>> jordan_normal_form(const Matrix<T> &A, const std::vector<std::pair<T, int>> &eigenvalues_list_typed) {
    if (A.empty()) {
        return {{}, {}};
    }
    size_t n = A.size();
    if (n > 0 && A[0].size() != (size_t)n) {
        throw std::invalid_argument("Input matrix must be square.");
    }

    Matrix<T> S_matrix(n, Vector<T>(n, T(0)));
    Matrix<T> J_matrix(n, Vector<T>(n, T(0)));
    std::vector<Vector<T>> S_columns_list;
    S_columns_list.reserve(n);
    int J_current_block_start_idx = 0;

    std::vector<std::pair<T, int>> sorted_eigenvalues = eigenvalues_list_typed;
    std::sort(sorted_eigenvalues.begin(), sorted_eigenvalues.end(),
              [](const auto &a, const auto &b) {
                  if constexpr (std::is_floating_point_v<T>) {
                      if (std::abs(a.first - b.first) > GLOBAL_EPSILON * 100) {
                          return a.first < b.first;
                      }
                  } else {
                      if (a.first != b.first) {
                          return a.first < b.first;
                      }
                  }
                  return a.second > b.second;
              });

    for (const auto &eig_pair : sorted_eigenvalues) {
        T lambda = eig_pair.first;
        int multiplicity = eig_pair.second;

        if (multiplicity == 0)
            continue;

        Matrix<T> B_operator = matrix_subtract(
            A, scalar_multiply_matrix(identity_matrix<T>(n), lambda));
        int p_max_chain_length =
            get_operator_index<T>(B_operator, multiplicity, n);

        if (p_max_chain_length == 0 && multiplicity > 0) {
            if (multiplicity > 0)
                p_max_chain_length = 1;
        }

        int vecs_found_for_this_lambda = 0;

        for (int k_chain_length = p_max_chain_length; k_chain_length >= 1;
             --k_chain_length) {
            if (vecs_found_for_this_lambda == multiplicity)
                break;

            Matrix<T> B_pow_k = matrix_power<T>(B_operator, k_chain_length);
            Matrix<T> B_pow_k_minus_1;
            if (k_chain_length == 1)
                B_pow_k_minus_1 = identity_matrix<T>(n);
            else
                B_pow_k_minus_1 =
                    matrix_power<T>(B_operator, k_chain_length - 1);

            std::vector<Vector<T>> basis_ker_Bk =
                get_nullspace_basis<T>(B_pow_k);
            if (basis_ker_Bk.empty() && k_chain_length > 0 &&
                multiplicity > vecs_found_for_this_lambda) {
                std::cerr
                    << "Warning: Ker(B^" << k_chain_length << ") for lambda "
                    << lambda
                    << " is empty or not found, but more vectors expected."
                    << std::endl;
                continue;
            }

            for (const auto &v_candidate_top : basis_ker_Bk) {
                if (is_zero_vector<T>(v_candidate_top))
                    continue;
                if (vecs_found_for_this_lambda == multiplicity)
                    break;

                Vector<T> first_vec_in_chain =
                    matrix_vector_multiply<T>(B_pow_k_minus_1, v_candidate_top);

                if (is_zero_vector<T>(first_vec_in_chain)) {
                    continue;
                }

                std::vector<Vector<T>> current_chain_vectors;
                current_chain_vectors.reserve(k_chain_length);
                for (int j = 0; j < k_chain_length; ++j) {
                    current_chain_vectors.push_back(matrix_vector_multiply<T>(
                        matrix_power<T>(B_operator, k_chain_length - 1 - j),
                        v_candidate_top));
                }

                if (are_vectors_linearly_independent<T>(
                        S_columns_list, current_chain_vectors)) {
                    S_columns_list.insert(S_columns_list.end(),
                                          current_chain_vectors.begin(),
                                          current_chain_vectors.end());
                    for (int i = 0; i < k_chain_length; ++i) {
                        J_matrix[J_current_block_start_idx + i]
                                [J_current_block_start_idx + i] = lambda;
                        if (i < k_chain_length - 1) {
                            J_matrix[J_current_block_start_idx + i]
                                    [J_current_block_start_idx + i + 1] = T(1);
                        }
                    }
                    J_current_block_start_idx += k_chain_length;
                    vecs_found_for_this_lambda += k_chain_length;
                }
            }
        }
        if (vecs_found_for_this_lambda != multiplicity) {
            std::cerr << "Warning: For eigenvalue " << lambda
                      << " (multiplicity " << multiplicity << "), "
                      << "expected " << multiplicity
                      << " basis vectors, but found "
                      << vecs_found_for_this_lambda << "." << std::endl;
        }
    }

    if (S_columns_list.size() != (size_t)n) {
        std::cerr << "Error: Expected to find " << n
                  << " basis vectors for S, but found " << S_columns_list.size()
                  << "." << std::endl;
    } else {
        std::cout << "Successfully found " << n
                  << " basis vectors for S matrix." << std::endl;
    }

    for (size_t j = 0; j < S_columns_list.size(); ++j) {
        if (S_columns_list[j].size() != (size_t)n) {
            throw std::runtime_error(
                "Internal error: S column vector has incorrect dimension "
                "during S_matrix assembly.");
        }
        for (int i = 0; i < n; ++i) {
            if (j < (size_t)n) {
                S_matrix[i][j] = S_columns_list[j][i];
            }
        }
    }
    return {J_matrix, S_matrix};
}



// A wrapper function for calling with a specific type
template <typename T>
void solve_for_type(int n, Matrix<T>& A_input_typed, std::vector<std::pair<T, int>>& eigenvalues_list_typed) {
    std::cout << "\n--- Computing Jordan Normal Form for type " << typeid(T).name() << " ---" << std::endl;
    try {
        auto [J_computed, S_computed] = jordan_normal_form<T>(A_input_typed, eigenvalues_list_typed);
        std::cout << "\n--- Results ---" << std::endl;
        print_matrix<T>("J_computed (Jordan Form)", J_computed);
        print_matrix<T>("S_computed (Transition Matrix)", S_computed);

        char verify_choice = 'n';
        std::cout << "Do you want to verify the result? (y/N): " << std::flush;
        std::cin >> verify_choice;
        if (verify_choice == 'y' || verify_choice == 'Y') {
            verify_jordan_decomposition<T>(A_input_typed, J_computed, S_computed);
        }

    } catch (const std::exception& e) {
        std::cerr << "An error occurred during computation: " << e.what() << std::endl;
    }
}



int main() {
    size_t n;
    std::cout << "Enter the size of the square matrix (n): " << std::flush;
    std::cin >> n;
    if (n == 0) {
        return 0;
    }

    int type_choice;
    std::cout << "Select data type for computations:\n";
    std::cout << "1. long double (for real numbers)\n";
    std::cout << "2. Modular arithmetic (over Z_p for prime p)\n";
    std::cout << "Enter choice (1 or 2): " << std::flush;
    std::cin >> type_choice;


    // ===== Long double branch: =====
    if (type_choice == 1) {
        std::cout << std::fixed << std::setprecision(DoublePrintPrecision);

        Matrix<long double> A(n, Vector<long double>(n));
        std::cout << "Enter input matrix A (long double):" << std::endl;
        for (size_t i = 0; i < n; ++i) {
            for(size_t j = 0; j< n; ++j) {
                std::cin >> A[i][j];
            }
        }

        std::vector<long double> raw_eigenvalues;
        std::cout << "Enter " << n << " eigenvalues (long double):" << std::endl;
        for(int i = 0; i < n; ++i) {
            long double val;
            std::cin >> val;
            raw_eigenvalues.push_back(val);
        }

        std::vector<std::pair<long double, int>> eigenvalues_to_multiplicity;
        std::sort(raw_eigenvalues.begin(), raw_eigenvalues.end());
        eigenvalues_to_multiplicity.push_back({raw_eigenvalues[0], 0});
        for (long double eig_val : raw_eigenvalues) {
            if (std::abs(eig_val - eigenvalues_to_multiplicity.back().first) < GLOBAL_EPSILON) {
                eigenvalues_to_multiplicity.back().second++;
            } else {
                eigenvalues_to_multiplicity.push_back({eig_val, 1});
            }
        }
        solve_for_type<long double>(n, A, eigenvalues_to_multiplicity);


    // ===== Modular arithmetic branch: =====
    } else if (type_choice == 2) {
        const int64_t P_MOD = PrimeForZpArithmetic;  // simple example of the prime module
        using MyModular = Modular<int64_t, P_MOD>;
        Matrix<MyModular> A(n, Vector<MyModular>(n));
        std::vector<std::pair<MyModular, int>> eigenvalues_to_multiplicity;

        std::cout << "Enter input matrix A (elements will be taken modulo " << P_MOD << "):" << std::endl;
        for (int i = 0; i < n; ++i) { 
            for(int j = 0; j < n; ++j) {
                std::cin >> A[i][j];
            }
        }

        std::vector<MyModular> raw_eigenvalues; 
        std::cout << "Enter " << n << " eigenvalues (will be taken modulo " << P_MOD << "):" << std::endl;
        for(int i = 0; i < n; ++i) {
            MyModular val;
            std::cin >> val;
            raw_eigenvalues.push_back(val);
        }

        std::sort(raw_eigenvalues.begin(), raw_eigenvalues.end()); 
        eigenvalues_to_multiplicity.push_back({raw_eigenvalues[0], 0});
        for (MyModular eig_val : raw_eigenvalues) {
            if (eig_val == eigenvalues_to_multiplicity.back().first) { 
                eigenvalues_to_multiplicity.back().second++;
            } else {
                eigenvalues_to_multiplicity.push_back({eig_val, 1});
            }
        }
        solve_for_type<MyModular>(n, A, eigenvalues_to_multiplicity);


    // ===== Incorrect choice: =====
    } else {
        std::cerr << "Invalid choice." << std::endl;
        return 1;
    }


    return 0;
}

