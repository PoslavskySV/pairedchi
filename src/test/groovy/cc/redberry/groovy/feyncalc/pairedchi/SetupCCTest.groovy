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

import cc.redberry.core.context.CC
import cc.redberry.core.context.OutputFormat
import cc.redberry.core.tensor.Product
import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.groovy.Redberry
import org.junit.Ignore
import org.junit.Test

import static cc.redberry.core.context.OutputFormat.SimpleRedberry
import static cc.redberry.core.context.OutputFormat.WolframMathematica
import static cc.redberry.core.tensor.Tensors.setAntiSymmetric
import static cc.redberry.core.utils.TensorUtils.info
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 03/04/15.
 */
class SetupCCTest {

    @Test
    public void testCalc_Sign() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC()
            def bottomSpin = 'scalar'
            def diags = stp.diagrams(bottomSpin) & (J = s)
            println diags.size()
            def g1 = 1, g2 = -1
            def results = []
            def pol = stp.setupPolarisations(g1, g2)

            for (def diag in diags) {
                def M2 = stp.calcProcess([diag], pol)
                M2 <<= 'g=1'.t &
                        's=900'.t &
                        'mc=15/10'.t &
                        'mb=45/10'.t &
                        't1=-684064/1000'.t &
                        't2=-110326/1000'.t &
                        'u1=-434806/100000'.t &
                        'u2=-499842/10000'.t
                results << M2
                println M2
            }
            stp.log "\n\n\n\n\n ******  $g1 $g2 done ****** \n\n\n\n\n\n\n"

            println results
        }
    }

    @Test
    public void testCalc() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC()
            def bottomSpin = 'scalar'
            def diags = stp.diagrams(bottomSpin)

            for (def g1 in [1, -1])
                for (def g2 in [1, -1]) {

                    def file = new File("/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/tmp-ss${g1}${g2}.redberry")
                    file.delete()

                    def pol = stp.setupPolarisations(g1, g2)
                    def M2 = stp.calcProcess(diags, pol)
                    M2 = M2 / stp.overallPolarizationFactor
                    file << M2.toString(OutputFormat.Redberry)

                    stp.log "\n\n\n\n\n ******************************  $g1 $g2 done ****************************** \n\n\n\n\n\n\n"
                }
        }
    }

    @Test
    public void testWardIdentities() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC()
            for (def bottomSpin in ['scalar', 'axial', 'tensor']) {

                def diags = stp.diagrams(bottomSpin)
                def diags_k1 = diags.collect { 'eps1_a[h1] = k1_a'.t >> it }
                def diags_k2 = diags.collect { 'eps2_a[h2] = k2_a'.t >> it }

                for (def g2 in [-1, 1]) {
                    def pol = stp.setupPolarisations(1, g2)
                    def M2 = stp.calcProcess(diags_k1, pol)
                    M2 <<= stp.mapleFactorTr
                    assert M2 == 0.t
                }

                for (def g1 in [-1, 1]) {
                    def pol = stp.setupPolarisations(g1, 1)
                    def M2 = stp.calcProcess(diags_k2, pol)
                    M2 <<= stp.mapleFactorTr
                    assert M2 == 0.t
                }

                println "\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~~ Done for $bottomSpin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n"
            }
        }
    }


    @Test
    public void testWardIdentities0() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC()
            def bottomSpin = 'scalar'
            def diags = stp.diagrams(bottomSpin).collect { ('eps1_a[h1] = k1_a'.t & 'eps2_a[h2] = k2_a'.t) >> it }
            def M2 = stp.calcProcess(diags)
            M2 <<= stp.mapleFactorTr
            assert M2 == 0.t
        }
    }

    @Test
    public void testName() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC()
            def bottomSpin = 'scalar'
            def diags = stp.diagrams(bottomSpin)

            for (def g1 in [1, -1]) {
                for (def g2 in [1, -1]) {
                    def file = new File("/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/ss${g1}${g2}.redberry")
                    file.delete()

                    def pol = stp.setupPolarisations(g1, g2)
                    def M2 = stp.calcProcess(diags, pol)
                    M2 = M2 / stp.overallPolarizationFactor
                    file << M2.toString(OutputFormat.Redberry)

                    stp.log "\n\n\n\n\n ******  $g1 $g2 done ****** \n\n\n\n\n\n\n"
                }
            }
        }
    }

    @Test
    public void testLevi() throws Exception {
        use(Redberry) {
            Setup stp = new Setup(false, true)
            def t = 'e_abcd * k1^a * k2^b * p1^c[charm] * p2^d[charm]'.t
            def t2 = t * t
            t2 <<= stp.leviSimplify & stp.fullSimplify & stp.massesSubs & stp.mFactor
            println t2
        }
    }

    @Test
    public void printAmpXYZ() throws Exception {

        use(Redberry) {
            CC.reset()
            CC.resetTensorNames(-2907357143612290431)

            SetupCC stp = new SetupCC()
            //axial: 21369 (not projected), 51048 (1, 1)
            def bottomSpin = 'scalar'

            def pol = Identity //stp.setupPolarisations(1, 1)
            def qXYZ = stp.setXYZ('bottom', 'z', 'z')
            def qpol = qXYZ['tr']


            def diags = stp.diagrams(bottomSpin)


            def i = 0
            def file = new File("/Users/poslavsky/Downloads/amps_inv_k1_z_z.txt")
            file.delete()
            file << "Gluons factor (inversed):"
            file << (stp.overallPolarizationFactor**(1.t / 2)).toString(WolframMathematica)
            file << "Bottom factor (inversed):"
            file << qXYZ['den'].toString(WolframMathematica)


            file << "\n\nMandelstam variables:\n\n"
            stp.mandelstam.each {
                file << "${it.toString(SimpleRedberry)};\n"
            }

            file << "\n\nSpinor structures:\n\n"
            stp.spinorStructures.each {
                file << "${it.toString(SimpleRedberry)};\n"
            }

            def tot = 0.t

            def stringify = 'p1_a[charm] = p1_a'.t.hold &
                    'p2_a[charm] = p2_a'.t.hold &
                    'eps1_a[h1] = eps1_a'.t.hold &
                    'eps2_a[h2] = eps2_a'.t.hold &
                    'p_a[bottom] = p_a'.t.hold &
                    'eps_{a}[h[bottom]] = epsP_a'.t.hold &
                    'eps_{ab}[h[bottom]] = epsP_ab'.t.hold




            def m2 = stp.calcProcess(diags.collect {('eps1_a[h1] = k1_a'.t & 'eps2_a[h2] = k2_a'.t) >> it }, qpol)

            println info(m2)
            println stp.mapleFactorTr >> m2

            return
            file << "\n\nAmplitudes:\n\n"
            for (def diag in diags) {
                stp.log "amp $i of ${diags.size()}"
                def amp = stp.calcAmplitude(diag, pol & qpol)
                amp <<= stringify


                def num = Numerator >> amp

//                def f = new File("/Users/poslavsky/Downloads/amp$i")
//                f.delete()
//                f << num.toString(OutputFormat.Redberry)

                num <<= ExpandTensors & EliminateMetrics & stp.mFactor

                amp = num / (Denominator >> amp)
                tot += amp

                amp <<= InvertIndices
                file << "AMP[$i][A, B] = ${amp.toString(WolframMathematica)};\n\n"
                ++i
            }

            println info(tot)

        }
    }


    @Test
    public void printAmp() throws Exception {

        use(Redberry) {
            CC.resetTensorNames(-2907357143612290431)
            SetupCC stp = new SetupCC()
            //axial: 21369 (not projected), 51048 (1, 1)
            def bottomSpin = 'axial'

            def diags = stp.diagrams(bottomSpin)
            def pol = stp.setupPolarisations(1, -1)

            def i = 0
            def file = new File("/Users/poslavsky/Downloads/amps_${bottomSpin}.txt")
            file.delete()

            file << "\n\nMandelstam variables:\n\n"
            stp.mandelstam.each {
                file << "${it.toString(SimpleRedberry)};\n"
            }

            file << "\n\nSpinor structures:\n\n"
            stp.spinorStructures.each {
                file << "${it.toString(SimpleRedberry)};\n"
            }

            def tot = 0.t

            def stringify = 'p1_a[charm] = p1_a'.t.hold &
                    'p2_a[charm] = p2_a'.t.hold &
                    'eps1_a[h1] = eps1_a'.t.hold &
                    'eps2_a[h2] = eps2_a'.t.hold &
                    'p_a[bottom] = p_a'.t.hold &
                    'eps_{a}[h[bottom]] = epsP_a'.t.hold &
                    'eps_{ab}[h[bottom]] = epsP_ab'.t.hold

            file << "\n\nAmplitudes:\n\n"
            for (def diag in diags) {
                stp.log "amp $i"
                def amp = stp.calcAmplitude(diag, pol)
                amp <<= stringify

                def num = Numerator >> amp
                num <<= ExpandTensors & EliminateMetrics & stp.mFactor

                amp = num / (Denominator >> amp)
                tot += amp

                amp <<= InvertIndices
                file << "AMP[$i][A, B] = ${amp.toString(WolframMathematica)};\n\n"
                ++i
            }

            println info(tot)
        }
    }

    @Test
    public void printMultTable() throws Exception {
        use(Redberry) {
            CC.resetTensorNames(-2907357143612290431)
            SetupCC stp = new SetupCC()
            def bottomSpin = 'scalar'

            def diags = stp.diagrams(bottomSpin)
            def pol = Identity //stp.setupPolarisations(1, 1)

            def i = 0
            def file = new File('/Users/poslavsky/Downloads/amps_axial.txt')
            file.delete()

            file << "\n\nMandelstam variables:\n\n"
            stp.mandelstam.each {
                file << "${it.toString(SimpleRedberry)};\n"
            }

            file << "\n\nSpinor structures:\n\n"
            stp.spinorStructures.each {
                file << "${it.toString(SimpleRedberry)};\n"
            }

            def tot = 0.t

            def stringify = 'p1_a[charm] = p1_a'.t.hold &
                    'p2_a[charm] = p2_a'.t.hold &
                    'eps1_a[h1] = eps1_a'.t.hold &
                    'eps2_a[h2] = eps2_a'.t.hold &
                    'p_a[bottom] = p_a'.t.hold &
                    'eps_{a}[h[bottom]] = epsP_a'.t.hold

            file << "\n\nAmplitudes:\n\n"
            for (def diag in diags) {
                stp.log "amp $i"
                def amp = stp.calcAmplitude(diag, pol)
                amp <<= stringify

                def num = Numerator >> amp
                num <<= ExpandTensors & EliminateMetrics & stp.mFactor
                num.each {
                    tot += it / (Denominator >> amp)
                }
            }

            tot <<= stp.mFactor
            println tot.size()

            def arr = []
            println info(tot)
        }
    }

    @Test
    public void testXXX() throws Exception {
        use(Redberry) {
            Setup stp = new Setup(false)
            stp.mandelstam.each {
                println it.toString(SimpleRedberry)
            }
        }

    }

    @Ignore
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
                str += rhs.toString(WolframMathematica)
                str += ';'
                println str
            }
        }
    }
}
