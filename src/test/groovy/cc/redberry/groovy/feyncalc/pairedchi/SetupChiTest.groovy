package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.core.context.CC
import cc.redberry.groovy.Redberry
import org.junit.Before
import org.junit.Test

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
}
