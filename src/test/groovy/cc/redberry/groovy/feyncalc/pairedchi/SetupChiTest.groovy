package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.core.context.CC
import cc.redberry.core.context.OutputFormat
import cc.redberry.groovy.Redberry
import org.junit.Before
import org.junit.Test

import static cc.redberry.groovy.RedberryStatic.ExpandAndEliminate
import static cc.redberry.groovy.RedberryStatic.Reduce

/**
 * Created by poslavsky on 02/04/15.
 */
class SetupChiTest {
    @Before
    public void before() throws Exception {
        CC.reset()
    }

    @Test
    public void testWardIdentities() throws Exception {
        use(Redberry) {
            SetupChi stp = new SetupChi()
            for (def charmSpin in ['scalar', 'axial', 'tensor'])
                for (def bottomSpin in ['scalar', 'axial', 'tensor']) {
                    def M = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
                    for (def i in [1, 2])
                        assert ({
                            def ward = ("eps${i}_i[h${i}] = k${i}_i".t
                                    & "k${3 - i}_i = p_i[bottom] - k${i}_i + p_i[charm]".t
                            ) >> M
                            ward <<= stp.fullSimplify & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 - s - t'.t & stp.mFactor
                            ward
                        }()) == 0.t
                }
        }
    }

    @Test
    public void test213() throws Exception {
        use(Redberry) {
            def x = 'x'.t
            def y = 'x + 1'.t
            def tr = 'x = 123'.t
            (x, y) = tr >> [x, y]
            println x
            println y
        }

    }

    @Test
    public void testPolarisations() throws Exception {
        use(Redberry) {
            SetupChi stp = new SetupChi()
            def eps1 = 'eps1_a = c1 * k1_a + c2 * k2_a + c3 * p_a[bottom]'.t
            def eps2 = 'eps2_a = c4 * e_abcd * k1^b * k2^c * p^d[bottom]'.t

            def epsPlus = 'epsPlus_a = (eps1_a + I * eps2_a)/2**(1/2)'.t
            def epsMinus = 'epsMinus_a = (eps1_a - I * eps2_a)/2**(1/2)'.t
            epsPlus <<= eps1 & eps2
            epsMinus <<= eps1 & eps2

            def eq = ['k1_a * eps1^a = 0'.t,
                      'k2_a * eps1^a = 0'.t,
                      'eps1_a * eps1^a = -1'.t,
                      'eps2_a * eps2^a = -1'.t,
                      'eps1_a * eps2^a = 0'.t
                      //'epsMinus_a * epsPlus_b + epsPlus_a * epsMinus_b = -g_ab'.t
            ]
            eq = (eps1 & eps2 & epsPlus & epsMinus & ExpandAndEliminate & stp.leviSimplify &
                    ExpandAndEliminate & stp.mandelstam & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s - t'.t) >> eq
            def options = [ExternalSolver: [
                    Solver: 'Mathematica',
                    Path  : '/Applications/Mathematica.app/Contents/MacOS']
            ]

            def s = Reduce(eq, ['c1', 'c2', 'c3', 'c4'].t, options)
            println s
            println s.size()
            s.each {
                println it
                println stp.wolframFactorTr >> it
            }

        }
    }

    @Test
    public void testPolarisations1() throws Exception {
        use(Redberry) {
            SetupChi stp = new SetupChi()
            def eps1 = 'eps1_a = s1 * p_a[charm] + s2 * p_a[bottom]'.t
            def eps2 = 'eps2_a = s3 * p_a[charm] + s4 * p_a[bottom] + s5 * k1_a'.t
            def eps0 = 'eps0_a = s6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]'.t

            def epsPlus = 'epsPlus_a = (eps1_a + I * eps2_a)/2**(1/2)'.t
            def epsMinus = 'epsMinus_a = (eps1_a - I * eps2_a)/2**(1/2)'.t
            epsPlus <<= eps1 & eps2
            epsMinus <<= eps1 & eps2

            def eq = ['p_a[charm] * eps1^a = 0'.t,
                      'p_a[charm] * eps2^a = 0'.t,
                      'eps1_a * eps2^a = 0'.t,
                      'eps1_a * eps1^a = -1'.t,
                      'eps2_a * eps2^a = -1'.t,
                      'eps0_a * eps0^a = -1'.t,
                      'p_a[charm] * eps0^a = 0'.t
                      //'epsMinus_a * epsPlus_b + epsPlus_a * epsMinus_b = -g_ab'.t
            ]
            eq = (eps1 & eps2 & eps0 & epsPlus & epsMinus & ExpandAndEliminate & stp.leviSimplify &
                    ExpandAndEliminate & stp.mandelstam & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s - t'.t) >> eq
            def options = [ExternalSolver: [
                    Solver: 'Mathematica',
                    Path  : '/Applications/Mathematica.app/Contents/MacOS']
            ]

            def s = Reduce(eq, ['s1', 's2', 's3', 's4', 's5', 's6'].t, options)
            println s.size()
            s[0].eachWithIndex { t, i ->
                i++
                println "def s$i = 's$i = ${stp.wolframFactorTr >> t}'.t"
            }

        }
    }


    @Test
    public void test231() throws Exception {
        use(Redberry) {
            SetupChi stp = new SetupChi()
            def c1 = 'c1 = -(t-4*mb**2)*(-t**2+4*t*mb**2-t*s-16*mc**2*mb**2+4*t*mc**2)**(-1/2)*s**(-1/2)'.t
            def c2 = 'c2 = (t-4*mc**2+s)*(-t**2+4*t*mb**2-t*s-16*mc**2*mb**2+4*t*mc**2)**(-1/2)*s**(-1/2)'.t
            def c3 = 'c3 = -(-t**2+4*t*mb**2-t*s-16*mc**2*mb**2+4*t*mc**2)**(-1/2)*s**(1/2)'.t
            def c4 = 'c4 = -2*(-s*(t**2-4*t*mb**2+t*s+16*mc**2*mb**2-4*t*mc**2))**(-1/2)'.t

            def eps1 = 'eps1_a = c1 * k1_a + c2 * k2_a + c3 * p_a[bottom]'.t
            def eps2 = 'eps2_a = c4 * e_abcd * k1^b * k2^c * p^d[bottom]'.t
            def epsPlus = 'epsPlus_a = (eps1_a + I * eps2_a)/2**(1/2)'.t
            def epsMinus = 'epsMinus_a = (eps1_a - I * eps2_a)/2**(1/2)'.t
            epsPlus <<= eps1 & eps2 & c1 & c2 & c3 & c4
            epsMinus <<= eps1 & eps2 & c1 & c2 & c3 & c4

            def sum = 'epsMinus_a * epsPlus_b + epsPlus_a * epsMinus_b = -g_ab'.t
            sum <<= epsPlus & epsMinus & ExpandAndEliminate & stp.leviSimplify & stp.mandelstam & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s -t'.t
            println stp.mFactor >> sum
            sum = 'epsMinus_a * epsMinus_b + epsPlus_a * epsPlus_b = -g_ab'.t
            sum <<= epsPlus & epsMinus & ExpandAndEliminate & stp.leviSimplify & stp.mandelstam & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s -t'.t
            println stp.mFactor >> sum
        }
    }

    @Test
    public void test234() throws Exception {
        use(Redberry) {
            SetupChi stp = new SetupChi()


            for (def charmSpin in ['axial', 'tensor'])
                for (def bottomSpin in ['axial', 'tensor']) {
                    def out = new FileOutputStream('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/' + charmSpin[0] + bottomSpin[0] + '.redberry')

                    def amps = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
                    def amp2 = stp.squareMatrixElement(amps)
                    def r = "charm${charmSpin}bottom${bottomSpin}Sum".t.eq(amp2)
                    out << r.toString(OutputFormat.Redberry)
                    out << '\n'

                    for (def g1 in [1, -1])
                        for (def g2 in [1, -1]) {
                            amps = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
                            amps <<= "h1 = $g1".t & "h2 = $g2".t
                            stp.log "Setting up gluon polarizations $g1 $g2 ..."
                            amps <<= stp.epsVals & stp.fullSimplifyE & stp.mFactor
                            stp.log '... done'
                            amp2 = stp.squareMatrixElement(amps)
                            amp2 <<= stp.mandelstam & stp.massesSubs
                            r = "charm${charmSpin}bottom${bottomSpin}[$g1, $g2]".t.eq(amp2)
                            out << r.toString(OutputFormat.Redberry)
                            out << '\n'
                        }
                }
        }
    }

}
