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

                    def diags = stp.diagrams(charmSpin, bottomSpin)
                    def diags_k1 = diags.collect { 'eps1_a[h1] = k1_a'.t >> it }
                    def diags_k2 = diags.collect { 'eps2_a[h2] = k2_a'.t >> it }

                    def M2 = '0'.t
                    for (def g2 in [-1, 1]) {
                        def pol = stp.setupPolarisations(1, g2)
                        M2 += stp.calcProcess(diags_k1, pol)
                    }
                    M2 <<= stp.mapleFactorTr
                    assert M2 == 0.t

                    M2 = '0'.t
                    for (def g1 in [-1, 1]) {
                        def pol = stp.setupPolarisations(g1, 1)
                        M2 += stp.calcProcess(diags_k2, pol)
                    }
                    M2 <<= stp.mapleFactorTr
                    assert M2 == 0.t
                }
        }
    }
}

