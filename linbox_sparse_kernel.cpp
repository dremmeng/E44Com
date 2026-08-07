// Compute a right nullspace basis of a sparse matrix over a prime field.
// Input:  prime rows cols nonzeros, followed by zero-based row col value triples.
// Output: rank nullity, followed by one basis vector per line.

#include <cstdint>
#include <iostream>
#include <stdexcept>

#include <givaro/modular.h>
#include <linbox/algorithms/gauss.h>
#include <linbox/matrix/dense-matrix.h>

int main() {
    try {
        std::int64_t prime;
        std::size_t rows, cols, nonzeros;
        if (!(std::cin >> prime >> rows >> cols >> nonzeros) || prime <= 1) {
            throw std::runtime_error("expected: prime rows cols nonzeros");
        }

        using Field = Givaro::Modular<std::int64_t>;
        Field field(prime);
        using Gauss = LinBox::GaussDomain<Field>;
        Gauss::Matrix input(field, rows, cols);

        for (std::size_t entry = 0; entry < nonzeros; ++entry) {
            std::size_t row, col;
            std::int64_t value;
            if (!(std::cin >> row >> col >> value) || row >= rows || col >= cols) {
                throw std::runtime_error("invalid sparse matrix entry");
            }
            Field::Element element;
            field.init(element, value);
            input.setEntry(row, col, element);
        }

        Gauss gauss(field);
        LinBox::DenseMatrix<Field> basis(field, cols, 0);
        gauss.nullspacebasisin(basis, input);

        const std::size_t nullity = basis.coldim();
        std::cout << (cols - nullity) << ' ' << nullity << '\n';
        for (std::size_t basis_col = 0; basis_col < nullity; ++basis_col) {
            for (std::size_t row = 0; row < cols; ++row) {
                if (row != 0) {
                    std::cout << ' ';
                }
                std::int64_t value;
                field.convert(value, basis.getEntry(row, basis_col));
                std::cout << value;
            }
            std::cout << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "linbox_sparse_kernel: " << error.what() << '\n';
        return 1;
    }
    return 0;
}