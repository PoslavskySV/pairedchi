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

/**
 * Created by poslavsky on 03/04/15.
 */
class ChiBCCTest {

    @Test
    public void test1123e4() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();
            def t = stp.get3GluonDiagrams('scalar')
            println TensorUtils.getAllDiffSimpleTensors(t)
            println t[0]
            println t[1]
            println t[2]
            t = stp.squareMatrixElement(t)
            t <<= stp.massesSubs
            println TensorUtils.getAllDiffSimpleTensors(t)
            println t.class
            println t[0]
        }
    }

    @Test
    public void test1123() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();
            for (def spin in ['scalar', 'axial', 'tensor']) {

                println 'GLUON'
                def t = stp.getGluonDiagrams(spin)
                println t.size()
                println t[0]
                println t[1]
                println t[2]

                println 'QUARK'
                t = stp.getQuarkDiagrams(spin)
                println t.size()
                println t[0]
                println t[1]
                println t[2]
                //todo assert content!
            }
        }
    }

    @Test
    public void testa121232312323() throws Exception {
        use(Redberry) {
            def spin = 'scalar'
            def path = '/Users/poslavsky/Projects/Mathematica/ChiPairedProduction/ChiBCC_ward_' + spin + '.m'
            def file = new File(path)
            if (file.exists()) {
                file.delete()
                file = new File(path)
            }

            ChiBCC stp = new ChiBCC();

            def diag1 = stp.getGluonDiagrams(spin)
            diag1 <<= 'eps2_a[h2] = k2_a'.t.hold & stp.mandelstam & stp.massesSubs

            def r1 = stp.massesSubs >> stp.squareMatrixElement(diag1)

            println 'exp'
            file << ('r1=' + r1.toString(OutputFormat.WolframMathematica) + ';')

//            r2 <<= ExpandAll
//            println TensorUtils.hasImaginaryPart(r2)
//
//            println 'r2'
//            //println r2
//
//
//            r2 <<= stp.mFactor
//            println r2
        }
    }

    @Test
    public void testWardForScalar() throws Exception {
        use(Redberry) {
            def spin = 'scalar'
            def path = '/Users/poslavsky/Projects/Mathematica/ChiPairedProduction/ChiBCC_ward_' + spin + '.m'
            def file = new File(path)
            if (file.exists()) {
                file.delete()
                file = new File(path)
            }

            ChiBCC stp = new ChiBCC();
            def methods = [0: stp.&getGluonDiagrams, 1: stp.&getQuarkDiagrams, 2: stp.&get3GluonDiagrams]
            def amps = new SumBuilder()

            for (int i = 0; i < 3; ++i) {
                def diag = methods[i](spin)
                diag <<= ('eps2_a[h2] = k2_a'.t.hold & stp.mandelstam & stp.massesSubs)
                amps << diag
            }
            def r = stp.mandelstam >> stp.squareMatrixElement(amps.build())
            file << "r = ${r.toString(OutputFormat.WolframMathematica)} ;"
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
