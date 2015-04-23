/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2015:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.core.context.OutputFormat
import cc.redberry.core.tensor.Product
import cc.redberry.core.tensor.SumBuilder
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.junit.Test

import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 03/04/15.
 */
class ChiBCCTest {

    @Test
    public void testGluonDiags_LorentzCharge() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();
            def qs = stp.getGluonDiagrams('scalar')
            recreateTempFile() << stp.squareMatrixElement(qs).toString(OutputFormat.WolframMathematica)
        }
    }

    @Test
    public void testQuarkDiagrams() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();
            def qs = stp.getQuarkDiagrams('scalar')
            assert !TensorUtils.getAllDiffSimpleTensors(qs).collect({ x -> x.name }).contains('d_ABC'.t.name)
        }
    }

    @Test
    public void testWardForScalar() throws Exception {
        use(Redberry) {
            def epss = 'eps2_a[h2] = k2_a'.t //& 'eps2_a[h2] = k2_a'.t

            def spin = 'scalar'
            def path = '/Users/poslavsky/Projects/Mathematica/ChiPairedProduction/ChiBCC_ward_1_' + spin + '.m'
            def file = new File(path)
            if (file.exists()) {
                file.delete()
                file = new File(path)
            }

            ChiBCC stp = new ChiBCC(epss);
            def methods = [1: stp.&getGluonDiagrams, 0: stp.&getQuarkDiagrams, 2: stp.&get3GluonDiagrams]
            def amps = new SumBuilder()

            for (int i = 0; i < 3; ++i) {
                def diag = methods[i](spin)

                diag <<= (epss & stp.mandelstam & stp.massesSubs)

                //recreateTempFile() << diag

                //def diag2 = (stp.mandelstam & stp.massesSubs) >> stp.squareMatrixElement(diag)
                //file << "r$i = ${diag2.toString(OutputFormat.WolframMathematica)} ;\n"

                amps << diag
            }
            def r = (stp.mandelstam & stp.massesSubs) >> stp.squareMatrixElement(amps.build())
            file << "r = ${r.toString(OutputFormat.WolframMathematica)} ;"
        }
    }

    @Test
    public void testWardForScalar_3quark() throws Exception {
        use(Redberry) {
            def epss = 'eps1_a[h1] = k1_a'.t & 'eps2_a[h2] = k2_a'.t
            def spin = 'scalar'
            ChiBCC stp = new ChiBCC(epss);

            def diag = stp.getQuarkDiagrams(spin)
            diag <<= (epss & stp.mandelstam & stp.massesSubs)

            def diag2 = stp.squareMatrixElement(diag)
            diag2 <<= stp.mandelstam & stp.massesSubs & stp.wolframFactorTr


            println diag2

            diag2 <<= 's = 123*x**2'.t & 't1 = 23*x**2'.t & 't2 = 13*x**2'.t & 'u1 = 12*x**2'.t & 'u2 = 143*x**2'.t & 'mc = 23*x'.t & 'mb = 43*x'.t
            diag2 <<= Factor
            println diag2
        }
    }

    @Test
    public void test123() throws Exception {
        use(Redberry) {
            def dims = getDimSubs()
            def epss = 'eps1_a[h1] = k1_a'.t & 'eps2_a[h2] = k2_a'.t
            ChiBCC stp = new ChiBCC(epss);

            println(dims >> stp.getQuarkDiagramsNotProjected())
            println(dims >> stp.getQuarkDiagrams('scalar'))
        }
    }

    @Test
    public void test123123() throws Exception {
        use(Redberry) {
            println getDimSubs() >> '(k1_a + k2_a)/(k1_a*n^a)'.t
        }
    }

    private static File getTempFile() {
        new File('/Users/poslavsky/Projects/Mathematica/ChiPairedProduction/temp')
    }

    private static File recreateTempFile() {
        File t = getTempFile()
        if (t.exists())
            t.delete()
        getTempFile()
    }

    private static getDimSubs() {
        use(Redberry) {

            Random rnd = new Random(213)
            def rd = { 0.1 + rnd.nextDouble() }
            def tr = Identity

            //all momentums
            ['k1_a', 'k1_a', 'k2_a', 'p_a[bottom]', 'p1_a[charm]', 'p2_a[charm]', 'q_a[bottom]', 'q_a[charm]'].each {
                tr &= "$it = ${rd()} * x * f_a".t
            }
            //all other vectors
            ['eps1_a[h]', 'eps2_a[h]', 'epsS_a[fl]'].each {
                tr &= "$it = ${rd()} * f_a".t
            }
            //all mandelstam
            ['s', 't1', 't2', 'u1', 'u2'].each {
                tr &= "$it = ${rd()} * x**2".t
            }
            //all masses
            ['mc', 'mb', 'm[charm]', 'm[bottom]'].each {
                tr &= "$it = ${rd()} * x".t
            }

            //matrices
            tr &= 'G_a = f_a'.t
            tr &= 'T_A = t_A'.t
            tr &= EliminateMetrics
            tr &= "v[p2_a[bottom]]*cu[p1_a[bottom]] = ${rd()}*x".t
            tr &= "v[p2_a[charm]]*cu[p1_a[charm]] = ${rd()}*x".t
            tr &= EliminateMetrics

            tr &= 'f_a*f^a = 1'.t
            return tr
        }
    }

    @Test
    public void test32() throws Exception {
        use(Redberry) {
            def epss = 'eps1_a[h1] := k1_a'.t & 'eps2_a[h2] := k2_a'.t
            ChiBCC stp = new ChiBCC();
            def expr = stp.getQuarkDiagrams('scalar')
            expr <<= epss & stp.mandelstam & stp.massesSubs
            //def expr = getTempFile().readLines()[0].t
            expr.each { println it }
//            expr <<= ExpandAll
//
//            println (stp.wolframFactorTr >> stp.squareMatrixElement(expr))
        }
    }

    @Test
    public void testBruteForce() throws Exception {
        use(Redberry) {

            for (def spin in ['scalar', 'axial', 'tensor']) {
                def path = '/Users/poslavsky/Projects/Mathematica/ChiPairedProduction/ChiBCC_' + spin + '.m'
                def file = new File(path)
                if (file.exists()) {
                    file.delete()
                    file = new File(path)
                }

                ChiBCC stp = new ChiBCC();
                def methods = [0: stp.&getGluonDiagrams, 1: stp.&getQuarkDiagrams, 2: stp.&get3GluonDiagrams]
                def amps = new SumBuilder()

                for (int i = 0; i < 3; ++i)
                    amps << methods[i](spin)

                def r = stp.massesSubs >> stp.squareMatrixElement(amps.build())
                file << "r = ${r.toString(OutputFormat.WolframMathematica)} ;"
            }
        }
    }

    @Test
    public void test12323() throws Exception {
        use(Redberry) {
            def spin = 'scalar'
            def path = '/Users/poslavsky/Projects/Mathematica/ChiPairedProduction/ChiBCC_' + spin + '.m'
            def file = new File(path)
            if (file.exists()) {
                file.delete()
                file = new File(path)
            }

            ChiBCC stp = new ChiBCC();
            //for (def spin in ['scalar', 'axial', 'tensor']) {


            def r1 = stp.squareMatrixElement(stp.getQuarkDiagrams(spin))
            def r2 = stp.squareMatrixElement(stp.getGluonDiagrams(spin))
            def r = r1 + r2

            println TensorUtils.hasImaginaryPart(r)
            println r.size()
            println 'r1'
            //r1 <<= stp.wolframFactorTr
            //println r1

            println 'r2'
            //r2 <<= stp.wolframFactorTr
            //println r2


            file << ('r1=' + r1.toString(OutputFormat.WolframMathematica) + ';')
            file << '\n'

            file << ('r2=' + r2.toString(OutputFormat.WolframMathematica) + ';')
            file << '\n'
        }
    }

    @Test
    public void printMandelstamForFeynCalc() throws Exception {
        use(Redberry) {
            Setup s = new Setup(false)
            s.mandelstam.transformations.each { expr ->
                Product lhs = expr[0]
                def rhs = s.massesSubs >> expr[1]
                def f = { x -> "Momentum[${x.stringName}]" }
                def str = ''
                str += "Pair[${f(lhs[0])},${f(lhs[1])}]"
                str += '='
                str += rhs.toString(OutputFormat.WolframMathematica)
                str += ';'
                println str
            }
        }
    }

    @Test
    public void testWardIdentities() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC()
            for (def charmSpin in ['scalar', 'axial', 'tensor'])
                for (def bottomSpin in ['scalar', 'axial', 'tensor']) {
                    def M = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
                    for (def i in [1, 2])
                        assert ({
                            def ward = ("eps${i}_i[h${i}] = k${i}_i".t
                                    & "k${3 - i}_i = p_i[bottom] - k${i}_i + p_i[charm]".t
                            ) >> M
                            ward <<= stp.fullSimplify & stp.massesSubs & stp.mFactor
                            ward
                        }()) == 0.t
                }
        }
    }
}
