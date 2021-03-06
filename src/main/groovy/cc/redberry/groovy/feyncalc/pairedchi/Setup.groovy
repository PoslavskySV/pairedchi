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

import cc.redberry.core.tensor.Expression
import cc.redberry.core.tensor.Product
import cc.redberry.core.tensor.SumBuilder
import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import com.wolfram.jlink.KernelLink
import com.wolfram.jlink.MathLinkFactory

import static cc.redberry.core.context.OutputFormat.WolframMathematica
import static cc.redberry.core.indices.IndexType.Matrix1
import static cc.redberry.core.indices.IndexType.Matrix2
import static cc.redberry.core.tensor.Tensors.setAntiSymmetric
import static cc.redberry.core.tensor.Tensors.setSymmetric
import static cc.redberry.groovy.RedberryPhysics.*
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 02/04/15.
 */
class Setup implements AutoCloseable {
    public static final Map spins = [scalar: 0, axial: 1, tensor: 2]

    /**
     * Project c cbar onto charmonium or not
     */
    final boolean projectCC
    final long startTime
    final boolean log

    /** Feynman rules and auxiliary tensors */
    public def J2, J4, PS, V, D, G, V3, FeynmanRules

    /** Spin projectors and polarisation sums */

    /**
     *  Map flavour -> substitution
     */
    public Map spinSingletProjector
    /**
     * Map total spin -> substitution
     */
    public Map totalSpinProjector
    /**
     * Sum over polarisations
     */
    public Transformation epsSum

    /** Transformations */
    public def mandelstam, massesSubs, momentums,
               dTrace, dTraceSimplify, uTrace, uSimplify, leviSimplify, fullSimplify, fullSimplifyE,
               conjugateSpinors, momentumConservation

    /**
     * Polarisations
     */
    public boolean calcPolarisations
    public def polarisations


    KernelLink mathematicaKernel
    public Transformation mFactor, wolframFactorTr, mSimplify, wolframSimplifyTr, wolframFactorSqrtTr

    public Setup(boolean projectCC) {
        this(projectCC, false)
    }

    public Setup(boolean projectCC, boolean log, boolean calcPolarisations = true) {
        this.startTime = System.currentTimeMillis()
        this.log = log
        this.projectCC = projectCC
        this.calcPolarisations = calcPolarisations
        init()
    }

    void init() {//basic initialization
        use(Redberry) {
            log 'Setting up basic tensors'

            //Setting up matrix objects
            defineMatrices 'T_A', Matrix2.matrix,//unitary matrices
                    'G_a', 'G5', Matrix1.matrix, //gamma matrices
                    'v[p_a]', 'u[p_a]', Matrix1.vector, Matrix2.vector, //final quark wave function
                    'cu[p_a]', 'cv[p_a]', Matrix1.covector, Matrix2.covector, //final antiquark wave function
                    'V_iA', Matrix1.matrix, Matrix2.matrix, //quark-gluon vertex
                    'D[p_m, mass]', Matrix1.matrix //quark propagator

            //Levi-Civita, SU(N) symmetric and structure tensors
            Quiet { setAntiSymmetric 'e_abcd' }
            Quiet { setSymmetric 'd_ABC' }
            Quiet { setAntiSymmetric 'f_ABC' }
            //polarization tensor
            Quiet { setSymmetric 'eps_ab[h]' }

            setupFeynRules()
            setupProjectors()
            setupTransformations()
            setupMathematica()
            if (calcPolarisations)
                setupPolarisations()

            log 'Setup finished'
        }
    }

    Tensor wolframFunc(String command, Map bindings) {
        use(Redberry) {
            bindings.each { k, v ->
                command = command.replace(k, (v.class == Expression ? v[1] : v).toString(WolframMathematica))
            }
            mathematicaKernel.evaluateToInputForm(command, 0).replace('^', '**').t
        }
    }

    void setupMathematica() {
        log 'Setting up Mathematica'
        String[] args;
        switch (os()) {
            case 'mac':
                args = ["-linkmode", "launch", "-linkname", "\"/Applications/Mathematica.app/Contents/MacOS/MathKernel\" -mathlink"];
                break
            case 'linux':
                args = ["-linkmode", "launch", "-linkname", "math -mathlink"]
                break
            default:
                throw new RuntimeException('No Windows')
        }

        mathematicaKernel = MathLinkFactory.createKernelLink(args)
        mathematicaKernel.discardAnswer();
        def mop = { tensor, func, bindings ->
            def t = wolframFunc(func, bindings)
            return tensor.class == Expression ? tensor[0].eq(t) : t
        }
        def wolframFactor = { tensor -> mop(tensor, 'Factor[expr]', [expr: tensor]) }
        def wolframSimplify = { tensor -> mop(tensor, 'Simplify[expr, Reals]', [expr: tensor]) }
        def wolframFactorSqrt = { tensor -> mop(tensor, '(expr)/.{Power[x_, y_] :> Power[Factor[Expand[x]], y]}', [expr: tensor]) }
        wolframFactorTr = wolframFactor as Transformation
        wolframSimplifyTr = wolframSimplify as Transformation
        wolframFactorSqrtTr = wolframFactorSqrt as Transformation
        mFactor = Factor[[FactorScalars: false, FactorizationEngine: wolframFactor]]
        mSimplify = Factor[[FactorScalars: false, FactorizationEngine: wolframSimplify]]
    }

    /**
     * Feynman rules
     */
    void setupFeynRules() {
        log 'Setting up Feynman rules'
        use(Redberry) {
            if (projectCC)
                PS = 'PS_mn[k_a] = -g_mn'.t
            else {
                PS = 'PS_mn[k_a] = -g_mn - (k_m*n_n + k_n*n_m)/(k_a*n^a) + n_a*n^a*k_m*k_n/(k_a*n^a)**2'.t
                //particular n
                PS <<= 'n_a = p_a[bottom]'.t
            }
            //Quark vertex
            V = 'V_mA = I*g*G_m*T_A'.t
            //Quark propagator
            D = 'D[p_m, mass] = I*(mass + p_m*G^m)/(p_m*p^m - mass**2)'.t
            //Gluon propagator
            G = PS >> 'G_mn[k_a] = I*PS_mn[k_a]/(k_a*k^a)'.t
            //3-gluon vertex
            V3 = 'V_{Aa Bb Cc}[k1_m, k2_m, k3_m] = g*f_{ABC}*(g_ab*(k1_c - k2_c) + g_bc*(k2_a - k3_a) + g_ca*(k3_b - k1_b))'.t
            //All Feynman rules in one transformation
            FeynmanRules = V & D & G & V3
        }
    }

    /**
     * Spin and polarisation sums
     */
    void setupProjectors() {
        log 'Setting up spin projectors and polarisation sums'
        use(Redberry) {
            // Auxiliary tensors
            J2 = 'J_ab[p_i, m] := p_a*p_b/m**2 - g_ab'.t
            J4 = 'J_abcd[p_i, m] := (J_ac[p_i, m]*J_bd[p_i, m] + J_ad[p_i, m]*J_bc[p_i, m])/2 - J_ab[p_i, m]*J_cd[p_i, m]/3'.t
            J4 <<= Expand

            /** quarkonia polarization projector */
            totalSpinProjector = [scalar: Identity, axial: Identity, tensor: Identity]
            /** sum over polarizations */

            epsSum = Identity
            conjugateSpinors = Identity
            spinSingletProjector = [:]
            momentums = [:]
            for (def fl in ['fl', 'charm', 'bottom']) { //for each flavour (`fl` denotes abstract flavour)
                //quarks momentums in terms of quarkonia momentum (q_i is relative momentum)
                momentums[fl] = "p2_m[$fl] = p_m[$fl]/2 + q_m[$fl]".t.hold & "p1_m[$fl] = p_m[$fl]/2 - q_m[$fl]".t.hold

                // set spin singlet projection
                spinSingletProjector[fl] = "v[p2_m[$fl]]*cu[p1_m[$fl]] = (p2_m[$fl]*G^m - m[$fl]) * epsS_m[$fl] * G^m * (p1_m[$fl]*G^m + m[$fl])".t

                // total spin projection
                totalSpinProjector['scalar'] &= "q_i[$fl]*epsS_j[$fl] = -J_ij[p_i[$fl], 2*m[$fl]]".t
                totalSpinProjector['axial'] &= "q_i[$fl]*epsS_j[$fl] = e_ijab * p^a[$fl] * eps^b[h[$fl]]".t
                totalSpinProjector['tensor'] &= "q_i[$fl]*epsS_j[$fl] = eps_ij[h[$fl]]".t

                //sum over polarizations for axial meson
                epsSum &= "eps_a[h[$fl]]*eps_b[h[$fl]] = J_ab[p_a[$fl], 2*m[$fl]]".t
                //sum over polarizations for tensor meson
                epsSum &= "eps_ab[h[$fl]]*eps_cd[h[$fl]] = J_abcd[p_a[$fl], 2*m[$fl]]".t

                // spinor sums
                epsSum &= "u[p1_m[$fl]]*cu[p1_m[$fl]] = m[$fl] + p1^m[$fl]*G_m".t
                epsSum &= "v[p2_m[$fl]]*cv[p2_m[$fl]] = -m[$fl] + p2^m[$fl]*G_m".t

                conjugateSpinors &= "v[p2_m[$fl]]*cu[p1_m[$fl]] = u[p1_m[$fl]]*cv[p2_m[$fl]]".t
            }

            //sum over polarizations for gluons
            epsSum &= PS >> 'eps1_a[h1]*eps1_b[h1] = PS_ab[k1_a]'.t
            epsSum &= PS >> 'eps2_a[h2]*eps2_b[h2] = PS_ab[k2_a]'.t
        }
    }

    /**
     * Setup all transformations
     */
    void setupTransformations() {
        log 'Setting up transformations'
        use(Redberry) {
            //Dirac, unitary traces and relative simplifications
            dTrace = DiracTrace['G_a']
            uTrace = UnitaryTrace['T_A', 'f_ABC', 'd_ABC', '3']
            uSimplify = UnitarySimplify['T_A', 'f_ABC', 'd_ABC', '3']
            leviSimplify = LeviCivitaSimplify.minkowski['e_abcd']

            //Mandelstam variables and mass shell
            if (projectCC)
                mandelstam = setMandelstam([k1_m: '0', k2_m: '0', 'p_m[charm]': '2*m[charm]', 'p_m[bottom]': '2*m[bottom]'])
            else
                mandelstam = setMandelstam5([k1_m: '0', k2_m: '0', 'p1_m[charm]': 'm[charm]', 'p2_m[charm]': 'm[charm]', 'p_m[bottom]': '2*m[bottom]'])

            /** Simplifications of polarizations */
            def simplifyPolarizations = 'eps1^a[h1] * k1_a = 0'.t &
                    'eps2^a[h2] * k2_a = 0'.t &
                    'eps^a[h[charm]] * p_a[charm] = 0'.t &
                    'eps^a[h[bottom]] * p_a[bottom] = 0'.t &
                    'eps^ab[h[charm]] * p_a[charm] = 0'.t &
                    'eps^ab[h[bottom]] * p_a[bottom] = 0'.t

            //substituting masses for factor
            massesSubs = 'm[charm] = mc'.t.hold & 'm[bottom] = mb'.t.hold

            /** Full simplification */
            def simplifyMetrics = EliminateMetrics & simplifyPolarizations & mandelstam &
                    'd^i_i = 4'.t & 'd^A_A = 8'.t & "d^i'_i' = 4".t & "d^A'_A' = 3".t
            if (projectCC)
                simplifyMetrics &= 'e_abcd * p^a[charm] * p^b[bottom] * k1^c * k2^d = 0'.t

            fullSimplify = simplifyMetrics &
                    ExpandAll[simplifyMetrics] & simplifyMetrics &
                    leviSimplify &
                    ExpandAll[simplifyMetrics] & simplifyMetrics

            fullSimplifyE = simplifyMetrics &
                    ExpandTensors[simplifyMetrics] & simplifyMetrics &
                    leviSimplify &
                    ExpandTensors[simplifyMetrics] & simplifyMetrics

            dTraceSimplify = DiracTrace[[Gamma: 'G_a', Simplifications: fullSimplify]]

            if (projectCC)
                momentumConservation = 'p_i[bottom] = k1_i + k2_i - p_i[charm]'.t.hold
            else
                momentumConservation = 'p_i[bottom] = k1_i + k2_i - p1_i[charm] - p2_i[charm]'.t.hold
        }
    }

    /**
     * First term in Taylor series
     * @param momentum differentiation variable
     */
    Transformation taylor(momentum) {
        use(Redberry) {
            { expr -> ((Differentiate[momentum.t] & momentum.t.eq(0.t)) >> expr) * momentum.t } as Transformation
        }
    }

    Map quarkoniaVertices = null

    /**
     * Calculating effective vertex g + g -> meson (scalar, axial, tensor)
     * @return map
     */
    Map effectiveQuarkoniaVertices() {
        if (quarkoniaVertices != null)
            return quarkoniaVertices;

        log 'Calculating effective vertex'
        use(Redberry) {
            def effVertices = [:]

            // main formula (two diagrams)
            def A = 'cu[p1_m[fl]]*(V_aA*D[p1_m[fl] - k1_m, m[fl]]*V_bB + V_bB*D[p1_m[fl] - k2_m, m[fl]]*V_aA)*v[p2_m[fl]]'.t
            // basic simplifications
            A <<= FeynmanRules & spinSingletProjector['fl'] &
                    dTrace & ExpandAndEliminate & uTrace & momentums['fl'] &
                    ExpandAll[EliminateMetrics] & EliminateMetrics &
                    'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t

            // Taylor expansion (scalar meson)
            def AScalar = ('q_i[fl] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[fl]'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & totalSpinProjector['scalar']
            ) >> A

            // Define
            effVertices['scalar'] = 'Ascalar_{aA bB}[fl, k1_m, k2_m]'.t.eq AScalar

            // Taylor expansion (axial meson)
            def AAxial = ('q_i[fl] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[fl]'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & totalSpinProjector['axial']
            ) >> A
            // Define
            effVertices['axial'] = 'Aaxial_{aA bB}[fl, k1_m, k2_m]'.t.eq AAxial

            // Taylor expansion (tensor meson)
            def ATensor = ('q_i[fl] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[fl]'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & totalSpinProjector['tensor']
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
            ) >> A
            // Define
            effVertices['tensor'] = 'Atensor_{aA bB}[fl, k1_m, k2_m]'.t.eq ATensor

            return quarkoniaVertices = Collections.unmodifiableMap(effVertices)
        }
    }

    Transformation pairVertex = null

    /**
     * Effective vertex g + g -> q + qbar
     */
    Transformation effectivePairVertex() {
        if (pairVertex != null)
            return pairVertex

        use(Redberry) {
            //two diagrams
            def B = 'cu[p1_m[fl]]*(V_aA*D[p1_m[fl] - k1_m, m[fl]]*V_bB + V_bB*D[p1_m[fl] - k2_m, m[fl]]*V_aA)*v[p2_m[fl]]'.t

            // Simplifying
            B <<= FeynmanRules & ExpandAll[EliminateMetrics] & EliminateMetrics &
                    'p1_m[fl]*p1^m[fl] = m[fl]**2'.t & 'p2_m[fl]*p2^m[fl] = m[fl]**2'.t

            return pairVertex = 'B_{aA bB}[fl, k1_m, k2_m]'.t.eq(B)
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// POLARIZATIONS ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    void setupPolarisations() {
        log 'Setting up polarizations'
        polarisations = Identity
        setupPolarizationCoefficients()
        setupGluonPolarisations()
        for (def fl in (projectCC ? ['charm', 'bottom'] : ['bottom']))
            setupQuarkoniaPolarisations(fl)
    }

    def polarizationCoefficients = Identity

    void setupPolarizationCoefficients() {
        use(Redberry) {
            if (true && projectCC) {
                //definitions for x#
                //def r1 = '16*mb**4 - 4*mb**2*s - 4*mb**2*t - 4*mb**2*u + t*u = x0**(-2)'.t
                //def r2 = '-4*mb**2 + 8*mb*mc - 4*mc**2 + s = x1**(-2)'.t
                //def r3 = '-4*mb**2 - 8*mb*mc - 4*mc**2 + s = x2**(-2)'.t
                //def r4 = '-128*mb**4*mc**2 - 128*mb**2*mc**4 + 16*mb**2*mc**2*s + 16*mb**4*t + 48*mb**2*mc**2*t - 4*mb**2*s*t - 4*mb**2*t**2 + 48*mb**2*mc**2*u + 16*mc**4*u - 4*mc**2*s*u - 4*mb**2*t*u - 4*mc**2*t*u + s*t*u - 4*mc**2*u**2 = x3**(-2)'.t
                //x0 -> x0/sqrt(s)
                //def x0 = 'x0**2 = s/(16*mb**4 - 4*mb**2*s - 4*mb**2*t - 4*mb**2*u + t*u)'.t
                def c1 = 'c1 = (-t+4*mb**2)*x0'.t
                def c2 = 'c2 = x0*(-u+4*mb**2)'.t
                def c3 = 'c3 = -x0*s'.t
                def c4 = 'c4 = -2*x0'.t
                def cs1 = 'cs1 = -(1/2)*mc**(-1)*x1*x2*(-s+4*mc**2+4*mb**2)'.t
                def cs2 = 'cs2 = -4*x1*mc*x2'.t
                def cs3 = 'cs3 = (-I)*(-1)**(-1/2)*x1*x2*x3*(48*mc**2*mb**2-4*mc**2*u-4*mb**2*u-8*t*mb**2+16*mb**4-4*s*mb**2+s*u)'.t
                def cs4 = 'cs4 = (-I)*(-1)**(-1/2)*x1*x2*(48*mc**2*mb**2-8*mc**2*u-4*t*mc**2+16*mc**4+t*s-4*t*mb**2-4*mc**2*s)*x3'.t
                def cs5 = 'cs5 = (-I)*x2**(-1)*x1**(-1)*(-1)**(-1/2)*x3'.t
                def cs6 = 'cs6 = -2*x3'.t
                def bs1 = 'bs1 = -4*mb*x1*x2'.t
                def bs2 = 'bs2 = -(1/2)*mb**(-1)*x1*x2*(-s+4*mc**2+4*mb**2)'.t
                def bs3 = 'bs3 = (-I)*(-1)**(-1/2)*x1*x2*x3*(48*mc**2*mb**2-4*mc**2*u-4*mb**2*u-8*t*mb**2+16*mb**4-4*s*mb**2+s*u)'.t
                def bs4 = 'bs4 = (-I)*(-1)**(-1/2)*x1*x2*(48*mc**2*mb**2-8*mc**2*u-4*t*mc**2+16*mc**4+t*s-4*t*mb**2-4*mc**2*s)*x3'.t
                def bs5 = 'bs5 = (-I)*x2**(-1)*x1**(-1)*(-1)**(-1/2)*x3'.t
                def bs6 = 'bs6 = -2*x3'.t
                polarizationCoefficients &= c1 & c2 & c3 & c4 & cs1 & cs2 & cs3 & cs4 & cs5 & cs6 & bs1 & bs2 & bs3 & bs4 & bs5 & bs6
                return
            }

            log 'solving equations for gluon polarizations'
            def eps1 = 'eps1_a = c1 * k1_a + c2 * k2_a + c3 * p_a[bottom]'.t
            def eps2 = 'eps2_a = c4 * e_abcd * k1^b * k2^c * p^d[bottom]'.t
            def eq = ["k1_a * eps1^a = 0".t,
                      "k2_a * eps1^a = 0".t,
                      'eps1_a * eps1^a = -1'.t,
                      'eps2_a * eps2^a = -1'.t,
                      'eps1_a * eps2^a = 0'.t]
            eq = (eps1 & eps2 & ExpandAndEliminate & leviSimplify &
                    ExpandAndEliminate & mandelstam & massesSubs) >> eq
            def options = [ExternalSolver: [
                    Solver: 'Mathematica',
                    Path  : '/Applications/Mathematica.app/Contents/MacOS']
            ]
            def solutions = Reduce(eq, ['c1', 'c2', 'c3', 'c4'].t, options)
            assert solutions.size() != 0
            polarizationCoefficients &= (wolframFactorTr & wolframFactorSqrtTr) >> solutions[0]

            for (def fl in (projectCC ? ['charm', 'bottom'] : ['bottom'])) {
                log "solving equations for ${fl}onium"
                def var = fl[0] + 's'
                eps1 = "eps1_a = ${var}1 * p_a[charm] + ${var}2 * p_a[bottom]".t
                eps2 = "eps2_a = ${var}3 * p_a[charm] + ${var}4 * p_a[bottom] + ${var}5 * k1_a".t
                def eps0 = "eps0_a = ${var}6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]".t

                if (!projectCC) {
                    def subs = 'p_a[charm] = p1_a[charm] + p2_a[charm]'.t.hold
                    eps1 <<= subs; eps2 <<= subs; eps0 <<= subs;
                }

                eq = ["p_a[$fl] * eps1^a = 0".t,
                      "p_a[$fl] * eps2^a = 0".t,
                      "p_a[$fl] * eps0^a = 0".t,
                      'eps1_a * eps2^a = 0'.t,
                      'eps1_a * eps1^a = -1'.t,
                      'eps2_a * eps2^a = -1'.t,
                      'eps0_a * eps0^a = -1'.t]
                eq = (eps1 & eps2 & eps0 & ExpandAndEliminate & leviSimplify &
                        ExpandAndEliminate & mandelstam & massesSubs) >> eq

                solutions = Reduce(eq, ["${var}1", "${var}2", "${var}3", "${var}4", "${var}5", "${var}6"].t, options)
                assert solutions.size() != 0
                polarizationCoefficients &= (wolframFactorTr & wolframFactorSqrtTr) >> solutions[0]
            }
        }
    }

    void setupGluonPolarisations() {
        log 'Setting up gluon polarisations'
        use(Redberry) {
            def eps1 = 'eps1_a = c1 * k1_a + c2 * k2_a + c3 * p_a[bottom]'.t
            def eps2 = 'eps2_a = c4 * e_abcd * k1^b * k2^c * p^d[bottom]'.t
            eps1 <<= polarizationCoefficients; eps2 <<= polarizationCoefficients
            for (def g in [1, 2]) {
                polarisations &= (eps1 & eps2) >> "eps${g}_a[1] = eps1_a".t//(eps1_a + I * eps2_a)/2**(1/2)".t
                polarisations &= (eps1 & eps2) >> "eps${g}_a[-1] = eps2_a".t//(eps1_a - I * eps2_a)/2**(1/2)".t
            }
        }
    }

    void setupQuarkoniaPolarisations(fl) {
        use(Redberry) {
            log "Setting up quarkonia  polarisations ($fl)"
            def var = fl[0] + 's'
            def eps1 = "eps1_a = ${var}1 * p_a[charm] + ${var}2 * p_a[bottom]".t
            def eps2 = "eps2_a = ${var}3 * p_a[charm] + ${var}4 * p_a[bottom] + ${var}5 * k1_a".t
            def eps0 = "eps0_a = ${var}6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]".t

            if (!projectCC) {
                def subs = 'p_a[charm] = p1_a[charm] + p2_a[charm]'.t.hold
                eps1 <<= subs; eps2 <<= subs; eps0 <<= subs;
            }

            //axial
            def epsPlus = "eps_a[$fl, 1] = (eps1_a + I * eps2_a)/2**(1/2)".t
            def epsMinus = "eps_a[$fl, -1] = (eps1_a - I * eps2_a)/2**(1/2)".t
            def epsZero = "eps_a[$fl, 0] = eps0_a".t

            for (def eps in [epsPlus, epsMinus, epsZero])
                polarisations &= (eps0 & eps1 & eps2 & polarizationCoefficients) >> eps

            def subs = Identity
            for (def sub in ['eps1_a * eps1_b',
                             'eps2_a * eps2_b',
                             'eps1_a * eps2_b',
                             'eps1_a * eps0_b',
                             'eps2_a * eps0_b',
                             'eps0_a * eps0_b'].t) {
                def rhs = sub
                rhs <<= eps0 & eps1 & eps2 & polarizationCoefficients & fullSimplifyE & massesSubs & mFactor
                subs &= sub.eq rhs
            }

            //tensor
            def tr = epsPlus & epsMinus & epsZero & ExpandAll & subs & ExpandTensors & mFactor
            polarisations &= tr >> "eps_ab[$fl, 2] = eps_a[$fl, 1] * eps_b[$fl, 1]".t
            polarisations &= tr >> "eps_ab[$fl, 1] = (eps_a[$fl, 1]*eps_b[$fl, 0] + eps_a[$fl, 0]*eps_b[$fl, 1])/2**(1/2)".t
            polarisations &= tr >> "eps_ab[$fl, 0] = 1/6**(1/2)*eps_a[$fl, 1]*eps_b[$fl, -1] + (2/3)**(1/2)*eps_a[$fl, 0]*eps_b[$fl, 0] + (1/6)**(1/2)*eps_a[$fl, -1]*eps_b[$fl, 1]".t
            polarisations &= tr >> "eps_ab[$fl, -1] = (eps_a[$fl, -1] * eps_b[$fl, 0] + eps_a[$fl, 0]*eps_b[$fl, -1])/2**(1/2)".t
            polarisations &= tr >> "eps_ab[$fl, -2] = eps_a[$fl, -1] * eps_b[$fl, -1]".t
        }
    }

    Transformation setPolarizations(def g1, def g2, def charmPol, def bottomPol,
                                    String bottomSpin, String charmSpin) {
        use(Redberry) {
            checkPol g1
            checkPol g2
            checkPolChi(charmPol, spins[charmSpin])
            checkPolChi(bottomPol, spins[bottomSpin])

            def subs = Identity
            if (g1 != null)
                subs &= "h1 = $g1".t
            if (g1 != null)
                subs &= "h2 = $g2".t
//            if (charmPol != null) {
//                subs &= "eps_a[h[charm]] = eps_a[charm, $charmPol]".t.hold
//                subs &= "eps_ab[h[charm]] = eps_ab[charm, $charmPol]".t.hold
//            }
//            if (bottomPol != null) {
//                subs &= "eps_a[h[bottom]] = eps_a[bottom, $bottomPol]".t.hold
//                subs &= "eps_ab[h[bottom]] = eps_ab[bottom, $bottomPol]".t.hold
//            }
//            if ([g1, g2, charmPol, bottomPol].any { it != null })
                subs &= polarisations

            return subs
        }
    }

    Tensor squareMatrixElement(Tensor matrixElement) {
        squareMatrixElement(matrixElement, null)
    }

    Tensor squareMatrixElement(Tensor matrixElement, String spins) {

        use(Redberry) {
            spins = spins == null ? '' : spins
            log('Squaring matrix element (size: ' + matrixElement.size() + ')')

            def invert = { expr ->
                (expr.indices.free.si % expr.indices.free.si.inverted) >> expr
            } as Transformation

            def conjugate = Conjugate & invert
            conjugate &= Reverse[Matrix1, Matrix2]
            conjugate &= conjugateSpinors

            //def cMatrixElement = conjugate >> matrixElement
            //return ParallelExpand.parallelExpand(matrixElement, cMatrixElement, this.&calcProduct as Transformation, { l -> this.log(l) })
            def totalTermsNumber = matrixElement.size() * matrixElement.size()
            def counter = 0

            def result = new SumBuilder()
            def percent = -1
            for (int i = 0; i < matrixElement.size(); ++i) {
                for (int j = 0; j < matrixElement.size(); ++j) {
                    def part = matrixElement[i] * (conjugate >> matrixElement[j])
                    if (part.class == Product)
                        part = calcProduct(part)

                    result << part

                    //print percentage
                    def p = (int) (100.0 * counter / totalTermsNumber)
                    if (p != percent) {
                        percent = p;
                        if (percent <= 10 || percent % 10 == 0)
                            log("$spins progress " + percent + "%")
                    }
                    ++counter
                }
            }
            return result.build()
        }
    }

    Tensor calcProduct(Tensor t) {
        use(Redberry) {
            if (t.class != Product)
                return (epsSum & uTrace & mandelstam & dTraceSimplify & fullSimplify & massesSubs & uSimplify) >> t
            Product part = t
            def indexless = part.indexlessSubProduct
            def tensor = part.dataSubProduct
            tensor <<= epsSum & uTrace & mandelstam & dTraceSimplify &
                    fullSimplify & massesSubs & uSimplify & massesSubs
            assert TensorUtils.isSymbolic(tensor)
            def res = indexless * tensor
            res <<= wolframFactorTr
            return res
        }
    }

    void log(String msg) {
        if (log) {
            int dur = (System.currentTimeMillis() - startTime) / 1000
            String elapsed = String.format("%d:%02d:%02d", (int) (dur / 3600), (int) (dur % 3600) / 60, (int) (dur % 60))
            println("[Redberry $elapsed] " + msg)
        }
    }

    @Override
    void close() throws Exception {
        if (mathematicaKernel != null)
            mathematicaKernel.close()
    }

    private static String os() {
        def os = System.getProperty('os.name').toLowerCase()
        if (os.contains('mac'))
            return 'mac'
        else if (os.contains('nux'))
            return 'linux'
        return null
    }

    static void checkPol(g) {
        if (g != null && g != 1 && g != -1)
            throw new IllegalArgumentException()
    }

    static void checkPolChi(g, j) {
        if (g != null && (g < -j || g > j))
            throw new IllegalArgumentException()
    }
}
