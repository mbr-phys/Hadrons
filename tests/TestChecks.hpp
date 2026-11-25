/*
 * TestChecks.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_tests_TestChecks_hpp_
#define Hadrons_tests_TestChecks_hpp_

#include <Hadrons/Global.hpp>

/*
 * Small helpers for tests that check module results numerically.
 *
 * Objects in the Hadrons environment are freed by the garbage collector
 * during and after Application::run(), so results cannot be fetched from
 * the environment afterwards. Instead every checked module writes its result
 * file (parameter `output` = stem(name)), and the test reads it back after
 * run() with the same ResultReader type (HDF5 or XML, as configured) and the
 * same top-level name the module passes to saveResult.
 *
 * Checks:
 *   exact   ||a - b|| <= tol max(||a||, ||b||)                 (counts as failure)
 *   vanish  ||a||     <= tol scale                              (counts as failure)
 *   report  relative difference printed only                    (approximate identities)
 * summary() prints the totals and returns EXIT_SUCCESS or EXIT_FAILURE.
 */

namespace HadronsTest
{
    typedef std::vector<Grid::Complex> CVec;

    // norm of a complex vector
    inline double norm(const CVec &a)
    {
        double n2 = 0.;

        for (auto &x: a)
        {
            n2 += std::norm(x);
        }

        return std::sqrt(n2);
    }

    // ||a - b|| / max(||a||, ||b||), 0 if both vanish, -1 if the sizes differ
    inline double relDiff(const CVec &a, const CVec &b)
    {
        if (a.size() != b.size())
        {
            return -1.;
        }

        CVec   d(a.size());
        double scale;

        for (unsigned int i = 0; i < a.size(); ++i)
        {
            d[i] = a[i] - b[i];
        }
        scale = std::max(norm(a), norm(b));

        return (scale == 0.) ? 0. : norm(d)/scale;
    }

    // y += c x (y is resized if empty)
    inline void axpy(CVec &y, const Grid::Complex c, const CVec &x)
    {
        if (y.empty())
        {
            y.assign(x.size(), Grid::Complex(0., 0.));
        }
        for (unsigned int i = 0; i < std::min(x.size(), y.size()); ++i)
        {
            y[i] += c*x[i];
        }
    }

    // flatten corr[momentum][time]
    inline CVec flat(const std::vector<CVec> &v)
    {
        CVec f;

        for (auto &x: v)
        {
            f.insert(f.end(), x.begin(), x.end());
        }

        return f;
    }

    class Checker
    {
    public:
        Checker(const std::string &dir, const unsigned int traj)
        : dir_(dir), traj_(traj)
        {}

        // output stem to pass as the module parameter `output`
        std::string stem(const std::string &name) const
        {
            return dir_ + "/" + name;
        }

        // file written by saveResult(stem(name), ...) for trajectory traj_
        // (same format as ModuleBase::resultFilename)
        std::string filename(const std::string &name) const
        {
            return stem(name) + "." + std::to_string(traj_) + "."
                   + Grid::Hadrons::resultFileExt;
        }

        // read back a result written by a module with top-level name tag
        template <typename T>
        T read(const std::string &name, const std::string &tag) const
        {
            T                                   res;
            Grid::Hadrons::ResultReader         reader(filename(name));

            Grid::read(reader, tag, res);

            return res;
        }

        void exact(const std::string &label, const CVec &a, const CVec &b,
                   const double tol)
        {
            const double r = relDiff(a, b);

            if (r < 0.)
            {
                fail(label, "size mismatch (" + std::to_string(a.size())
                     + " vs " + std::to_string(b.size()) + ")");
            }
            else if (a.empty())
            {
                fail(label, "empty data");
            }
            else
            {
                record(label, r, tol, (r <= tol), "rel. diff.");
            }
        }

        void vanish(const std::string &label, const CVec &a, const double scale,
                    const double tol)
        {
            const double r = (scale > 0.) ? norm(a)/scale : norm(a);

            if (a.empty())
            {
                fail(label, "empty data");
            }
            else
            {
                record(label, r, tol, (r <= tol), "rel. norm");
            }
        }

        void report(const std::string &label, const CVec &a, const CVec &b)
        {
            const double r = relDiff(a, b);

            std::cout << "[INFO] " << label << ": approximate identity, rel. diff. = "
                      << r << " (reported only)" << std::endl;
        }

        void fail(const std::string &label, const std::string &why)
        {
            nFail_++;
            std::cout << "[FAIL] " << label << ": " << why << std::endl;
        }

        int summary(void) const
        {
            std::cout << "==== " << nPass_ << " check(s) passed, " << nFail_
                      << " failed ====" << std::endl;

            return (nFail_ == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
        }
    private:
        void record(const std::string &label, const double r, const double tol,
                    const bool ok, const std::string &what)
        {
            if (ok)
            {
                nPass_++;
            }
            else
            {
                nFail_++;
            }
            std::cout << (ok ? "[PASS] " : "[FAIL] ") << label << ": " << what
                      << " = " << r << " (tol " << tol << ")" << std::endl;
        }
    private:
        std::string  dir_;
        unsigned int traj_, nPass_{0}, nFail_{0};
    };
}

#endif // Hadrons_tests_TestChecks_hpp_
