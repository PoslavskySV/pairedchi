package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.core.context.CC
import cc.redberry.core.context.OutputFormat
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.junit.Before
import org.junit.Test

import static cc.redberry.groovy.RedberryStatic.*

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
    public void test123() throws Exception {
        use(Redberry) {
            SetupChi.calc(1, 1, 2, 2, 'tensor', 'tensor', new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/tt1122.redberry'))
        }
    }

    @Test
    public void test11() throws Exception {
        SetupChi stp = new SetupChi()
        def pols = stp.setPolarizations(1, 1, 2, 2, 'tensor', 'tensor')
        def amps = stp.setupFeynmanDiagrams('tensor', 'tensor', pols)
        println TensorUtils.getAllDiffSimpleTensors(amps)
        println amps
    }

    @Test
    public void testSingleDiag() throws Exception {
        use(Redberry) {
            SetupChi stp = new SetupChi()
            def charmSpin = 'tensor'
            def bottomSpin = 'tensor'
            def pols = stp.setPolarizations(1, 1, 0, 0, 'tensor', 'tensor')

            def Ma = "eps1^a[h1] * A${charmSpin}_{aA cC}[charm, k1_i, -k1_i + p_i[charm]] * G^cd[k1_i - p_i[charm]] * g^CD * eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]".t
            def Mb = "eps1^a[h1] * A${bottomSpin}_{aA cC}[bottom, k1_i, -k1_i + p_i[bottom]] * G^cd[k1_i - p_i[bottom]] * g^CD * eps2^b[h2] * A${charmSpin}_{dD bB}[charm, -k2_i + p_i[charm], k2_i]".t

            def Mt = 0.t
            for (def M in [Ma, Mb]) {
                M <<= stp.FeynmanRules & stp.qVertices & stp.mandelstam
                def n = Numerator >> M
                n <<= stp.fullSimplify & stp.massesSubs & stp.mFactor
                n <<= pols & stp.fullSimplifyE & stp.massesSubs & stp.mFactor
                def d = Denominator >> M
                d <<= stp.fullSimplify & stp.massesSubs & stp.mFactor
                Mt += n / d
            }

            Mt <<= 'x0 = 1'.t
            def Mtc = (Conjugate & InvertIndices) >> Mt
            def x02 = 'x0**2 = s/(16*mb**4 - 4*mb**2*s - 4*mb**2*t - 4*mb**2*u + t*u)'.t
            def res = 0.t
            for (int i = 0; i < Mt.size(); ++i)
                for (int j = 0; j < Mt.size(); ++j) {
                    println(i + ' ' + j)
                    def term = Mt[i] * Mtc[j]

                    def d = Denominator >> term
                    println d
                    def n = Numerator >> term
                    n <<= ExpandTensors[EliminateMetrics] & stp.epsSum & stp.mandelstam &
                            stp.fullSimplifyE & stp.massesSubs
                    //d <<= stp.fullSimplify & stp.massesSubs & stp.mFactor
                    res += n / d
                }

            res <<= 'u = 4*mc**2 + 4*mb**2 - s - t'.t & stp.wolframFactorTr
            println TensorUtils.info(res)
            println res.toString(OutputFormat.WolframMathematica)
        }
    }

}

