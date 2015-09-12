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
import cc.redberry.core.utils.THashMap
import cc.redberry.groovy.Redberry
import org.junit.Ignore
import org.junit.Test

import static cc.redberry.core.context.OutputFormat.SimpleRedberry
import static cc.redberry.core.context.OutputFormat.WolframMathematica
import static cc.redberry.core.utils.TensorUtils.*
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 03/04/15.
 */
class SetupCCTest {

    public static Tensor replaceTensors(Map map, Tensor expr0) {
        use(Redberry) {
            def expr = (map.collect { "$it.key = $it.value".t }) >>> expr0

            def tts = []
            expr.parentAfterChild { t ->
                if (t.class == Product && t.indices.size() != 0)
                    tts << t.dataSubProduct
            }

            def var = map.size()
            for (def ts in tts) {
                def subs = "$ts = var${var + 1}".t
                def mod = subs >>> expr
                if (mod != expr) {
                    map[ts] = "var${++var}".t
                    expr = mod
                }
            }
            if (expr == expr0)
                return expr
            else return replaceTensors(map, expr)
        }
    }

    private static boolean freeQ(Tensor t, Tensor m) {
        use(Redberry) {
            t.parentAfterChild {
                if ((it % m).exists || (m % it).exists)
                    return false
            }
            return true
        }
    }

    @Test
    public void testCalc_Structs() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC(true)
            def bottomSpin = 'tensor'
            def diags = stp.diagrams(bottomSpin)
            def qpol = stp.setXYZ('bottom', 'x', 'x')['tr']
            def pol = stp.setupPolarisations(1, -1)

            def suntr = Identity
            suntr &= 'L1_AB = L1*g_AB'.t
            suntr &= 'L2_AB = L1*tt_AB'.t
            suntr &= 'L3_AB = L3*g_AB'.t
            suntr &= 'L4_AB = L3*tt_AB'.t
            suntr &= 'L5^i_AB = L5^i*g_AB'.t
            suntr &= 'L6^i_AB = L5^i*tt_AB'.t
            suntr &= 'L7^i_AB = L7^i*g_AB'.t
            suntr &= 'L8^i_AB = L7^i*tt_AB'.t
            suntr &= 'L9^ij_AB = L9^ij*g_AB'.t
            suntr &= 'L10^ij_AB = L9^ij*tt_AB'.t
            suntr &= 'L11^ij_AB = L11^ij*g_AB'.t
            suntr &= 'L12^ij_AB = L11^ij*tt_AB'.t
            suntr &= 'L13^ijk_AB = L13^ijk*g_AB'.t
            suntr &= 'L14^ijk_AB = L13^ijk*tt_AB'.t
            suntr &= 'L15^ijk_AB = L15^ijk*g_AB'.t
            suntr &= 'L16^ijk_AB = L15^ijk*tt_AB'.t


            def tot = 0.t
            def map = [:]
            for (def diag in diags) {
                def amp = stp.calcAmplitude(diag, qpol & pol)
                def num = Numerator >> amp, den = Denominator >> amp

                num = suntr >>> num
                num <<= ExpandTensorsAndEliminate
                num <<= 'e_abcd*k1^a*k2^b*p1^c[charm]*p2^d[charm] = lc'.t

                num *= 'f^AB'.t
                num <<= ExpandTensorsAndEliminate & stp.simplifyMetrics
                num <<= 'f^A_A = SUNN'.t
                num <<= 'f^AB*tt_AB = SUNT1'.t
                num <<= 'f^AB*tt_BA = SUNT2'.t
                num <<= ExpandTensorsAndEliminate
                num <<= stp.wFactor

                tot += replaceTensors(map, num) / den

                if (!isSymbolic(tot)) {
                    tot.parentAfterChild {
                        if (it.class == Product && it.indices.size() != 0) {
                            println 'beda \n\n'
                            println it
                            println '---beda \n\n'
                        }
                    }
                }
            }

            println map.size()
            println '\n\n\n'
            map.each { k, v ->
                println "$k   =   $v"
            }
            println '\n\n\n'
            println info(tot)
        }
    }

    @Test
    public void testCalc_Sign() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC()
            def bottomSpin = 'scalar'
            def diags = stp.diagrams(bottomSpin)
            println diags.size()
            def g1 = 1, g2 = -1
            def results = []
            def pol = stp.setupPolarisations(g1, g2)

            for (def diag in diags) {
                def M2 = stp.calcProcess([diag], pol)
                M2 <<= 'g=1'.t &
                        's=900*x**2'.t &
                        'mc=15/10*x'.t &
                        'mb=45/10*x'.t &
                        't1=-684064/1000*x**2'.t &
                        't2=-110326/1000*x**2'.t &
                        'u1=-434806/100000*x**2'.t &
                        'u2=-499842/10000*x**2'.t
                results << M2
//                assert isRealPositiveNumber(M2)
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
                for (def g2 in [-1, 1]) {

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
    public void testXYZWard() throws Exception {
        use(Redberry) {
            SetupCC stp = new SetupCC(true)
            def qXYZ = stp.setXYZ('bottom', 'x', 'y')
            def pol = stp.setupPolarisations(-1, 1)
            def qpol = qXYZ['tr']

            def diags = stp.diagrams('scalar')
//            diags = diags.collect { ('eps1_a[h1] = k1_a'.t) >> it }

            def m2 = stp.calcProcess(diags, pol & qpol)
            println info(m2)

            m2 <<= stp.mapleFactorTr

            println m2
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




            def m2 = stp.calcProcess(diags.collect { ('eps1_a[h1] = k1_a'.t & 'eps2_a[h2] = k2_a'.t) >> it }, qpol)

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

    @Test
    public void testX() throws Exception {
        use(Redberry) {
            def path = '/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/'
            def map = new THashMap()
            for (def i in 0..35) {
                def f = new File(path + "xyz__$i")

                def st = false
                for (def line in f.readLines()) {
                    if (st && (line.startsWith('Diagrams') || line.isEmpty()))
                        break
                    else if (line.startsWith('Lorentz structures'))
                        st = true
                    else if (st) {
                        def e = line.t
                        map[e[1]] = e[0]
                    }
                }
            }

            println map.size()

            map.each { k, v ->
                println "$k -> $v"
            }
        }
    }

    @Test
    public void test333() throws Exception {
        use(Redberry) {
            def d0 = '(4*mb**2-u2+2*mc**2-s-u1)**(-2)*(var2*((64*I)*TBA*g**4*s**2*(t2*u1+s*t1-4*mc**2*s+s*u1+4*mc**4-2*mc**2*t1-2*mc**2*u1+s**2-2*u2*mc**2-4*s*mb**2+t1*u1+t2*u2-2*t2*mc**2+u2*s+u2*t1+t2*s)**2*mb**2*(mc**2-t1)*(-t2*s**2+t2*u2*t1-3*u2*mc**2*t1+mc**2*t1*u1+s*t1*u1+s**2*t1+s*t1**2+2*mc**4*t1+4*t2*mc**2*s-mc**2*t1**2-4*s*t1*mb**2-2*mc**4*u1-t2*s*u1+2*u2*mc**4+3*t2*mc**2*u1-t2*u2*mc**2-2*t2*mc**4-u1*t2**2-s*t2**2+4*t2*s*mb**2-t2*t1*u1-4*mc**2*s*t1+mc**2*t2**2+u2*s*t1-t2*u2*s+u2*t1**2)+(-64*I)*g**4*s**2*(-t2+mc**2)*TAB*(t2*u1+s*t1-4*mc**2*s+s*u1+4*mc**4-2*mc**2*t1-2*mc**2*u1+s**2-2*u2*mc**2-4*s*mb**2+t1*u1+t2*u2-2*t2*mc**2+u2*s+u2*t1+t2*s)**2*mb**2*(t2*s**2-t2*u2*t1+3*u2*mc**2*t1-mc**2*t1*u1-s*t1*u1-s**2*t1-s*t1**2-2*mc**4*t1-4*t2*mc**2*s+mc**2*t1**2+4*s*t1*mb**2+2*mc**4*u1+t2*s*u1-2*u2*mc**4-3*t2*mc**2*u1+t2*u2*mc**2+2*t2*mc**4+u1*t2**2+s*t2**2-4*t2*s*mb**2+t2*t1*u1+4*mc**2*s*t1-mc**2*t2**2-u2*s*t1+t2*u2*s-u2*t1**2))+var1*((64*I)*(mc**2*u2**2+t2*s**2+3*u2*mc**2*t1+mc**2*t1*u1+t2*u2*u1-4*s*u1*mb**2+t2*u1**2+4*u2*s*mb**2+2*s*t1*u1-u2**2*s+s**2*t1+4*mc**4*s-u2*t1*u1+2*s**2*u1+s*u1**2-4*mc**2*s**2-u2**2*t1-2*mc**4*t1-2*t2*mc**2*s+2*mc**4*u1+s**3+2*t2*s*u1-2*u2*mc**4-mc**2*u1**2-4*s**2*mb**2-3*t2*mc**2*u1-t2*u2*mc**2+2*t2*mc**4-2*mc**2*s*t1-6*mc**2*s*u1+2*u2*mc**2*s)*g**4*s**2*(-t2+mc**2)*TAB*(t2*u1+s*t1-4*mc**2*s+s*u1+4*mc**4-2*mc**2*t1-2*mc**2*u1+s**2-2*u2*mc**2-4*s*mb**2+t1*u1+t2*u2-2*t2*mc**2+u2*s+u2*t1+t2*s)**2*mb**2+(-64*I)*TBA*g**4*s**2*(-mc**2*u2**2+t2*s**2-3*u2*mc**2*t1-mc**2*t1*u1-t2*u2*u1+4*s*u1*mb**2-t2*u1**2-4*u2*s*mb**2+u2**2*s+s**2*t1+4*mc**4*s+u2*t1*u1-s*u1**2-4*mc**2*s**2+u2**2*t1+2*mc**4*t1-2*t2*mc**2*s-2*mc**4*u1+s**3+2*u2*mc**4+mc**2*u1**2-4*s**2*mb**2+3*t2*mc**2*u1+t2*u2*mc**2-2*t2*mc**4-2*mc**2*s*t1+2*mc**2*s*u1+2*u2*s**2+2*u2*s*t1+2*t2*u2*s-6*u2*mc**2*s)*(t2*u1+s*t1-4*mc**2*s+s*u1+4*mc**4-2*mc**2*t1-2*mc**2*u1+s**2-2*u2*mc**2-4*s*mb**2+t1*u1+t2*u2-2*t2*mc**2+u2*s+u2*t1+t2*s)**2*mb**2*(mc**2-t1))+(128*I)*L1*g**4*s**3*mc*(-t2+mc**2)*TAB*(t2*u1+s*t1-4*mc**2*s+s*u1+4*mc**4-2*mc**2*t1-2*mc**2*u1+s**2-2*u2*mc**2-4*s*mb**2+t1*u1+t2*u2-2*t2*mc**2+u2*s+u2*t1+t2*s)**3*mb**2+(128*I)*TBA*L1*g**4*s**3*mc*(t2*u1+s*t1-4*mc**2*s+s*u1+4*mc**4-2*mc**2*t1-2*mc**2*u1+s**2-2*u2*mc**2-4*s*mb**2+t1*u1+t2*u2-2*t2*mc**2+u2*s+u2*t1+t2*s)**3*mb**2*(mc**2-t1))*(-t2+mc**2)**(-2)*(mc**2-t1)**(-2)*(-u2+2*mc**2-s-u1)**(-4)'.t
            def d2 = '((-128*I)*TBA*g**4*s**3*(-t2*u1-s*t1+4*mc**2*s-s*u1-4*mc**4+2*mc**2*t1+2*mc**2*u1-s**2+2*u2*mc**2+4*s*mb**2-t1*u1-t2*u2+2*t2*mc**2-u2*s-u2*t1-t2*s)**3*mb**2+(128*I)*g**4*s**3*TAB*(-t2*u1-s*t1+4*mc**2*s-s*u1-4*mc**4+2*mc**2*t1+2*mc**2*u1-s**2+2*u2*mc**2+4*s*mb**2-t1*u1-t2*u2+2*t2*mc**2-u2*s-u2*t1-t2*s)**3*mb**2)*(4*mb**2-u2+2*mc**2-s-u1)**(-2)*(-t2+4*mb**2-u2+4*mc**2-s-t1-u1)**(-2)*var1*(-u2+2*mc**2-s-u1)**(-4)'.t

            def subs = 'TBA=1'.t &
                    'TAB=2'.t &
                    'GAB=1'.t &
                    'g=1'.t &
                    's=900*x**2'.t &
                    'mc=15/10*x'.t &
                    'mb=45/10*x'.t &
                    't1=-684064/1000*x**2'.t &
                    't2=-110326/1000*x**2'.t &
                    'u1=-434806/100000*x**2'.t &
                    'u2=-499842/10000*x**2'.t & ExpandAll

            println subs >> d0
            println subs >> d2
        }
    }
}
