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
package cc.redberry.groovy.feyncalc.pairedchi;

import cc.redberry.core.context.CC;
import cc.redberry.core.indices.IndicesFactory;
import cc.redberry.core.tensor.ApplyIndexMapping;
import cc.redberry.core.tensor.Sum;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.Tensors;
import cc.redberry.core.tensor.random.RandomTensor;
import cc.redberry.core.transformations.*;
import cc.redberry.core.transformations.expand.ExpandTransformation;
import cc.redberry.core.transformations.expand.ExpandUtils;
import cc.redberry.core.utils.TensorUtils;
import cc.redberry.physics.feyncalc.LeviCivitaSimplifyTransformation;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;

import static cc.redberry.core.tensor.Tensors.*;

/**
 * Created by poslavsky on 24/04/15.
 */
@Ignore
public class ParallelExpandTest {
    @Test
    public void test1() throws Exception {
        CC.resetTensorNames(123);

        setAntiSymmetric("e_abcd");
        setAntiSymmetric("D_ac");
        setAntiSymmetric("B_abc");
        setAntiSymmetric("A_abc");

        RandomTensor rnd = new RandomTensor(false);
        rnd.reset(123);

        rnd.addToNamespace(parse("F_a"));
        rnd.addToNamespace(parse("A_ab"));
        rnd.addToNamespace(parse("B_abc"));
        rnd.addToNamespace(parse("D_ac"));
        rnd.addToNamespace(parse("g_ac"));
        rnd.addToNamespace(parse("e_abcd"));


        Tensor t1 = rnd.nextSum(80, 8, IndicesFactory.EMPTY_INDICES);
        Tensor t2 = rnd.nextSum(80, 8, IndicesFactory.EMPTY_INDICES);
        Transformation tr = new TransformationCollection(
                EliminateMetricsTransformation.ELIMINATE_METRICS,
                Tensors.parseExpression("A_ab*B^bac = T^c"),
                Tensors.parseExpression("A_ab*A^ba = xx"),
                Tensors.parseExpression("D_ab*D^ba = yy"),
                EliminateDueSymmetriesTransformation.ELIMINATE_DUE_SYMMETRIES,
                new LeviCivitaSimplifyTransformation(parseSimple("e_abcd"), true),
                ExpandAndEliminateTransformation.EXPAND_AND_ELIMINATE,
                Tensors.parseExpression("A_ab*B^bac = T^c"),
                Tensors.parseExpression("A_ab*A^ba = xx"),
                Tensors.parseExpression("D_ab*D^ba = yy")
        );

        System.out.println("Cold JVM: ");

        List<Tensor> results = new ArrayList<>();
        long start = System.currentTimeMillis();
        results.add(Transformation.Util.applyUntilUnchanged(new ExpandTransformation(tr).transform(multiplyAndRenameConflictingDummies(t1, t2)), tr));
        System.out.println("Normal: " + (System.currentTimeMillis() - start));

        for (int i = 1; i < 5; ++i) {
            start = System.currentTimeMillis();
            ParallelExpand pe = new ParallelExpand((Sum) t1, (Sum) t2, i, tr, null);
            results.add(Transformation.Util.applyUntilUnchanged(pe.expand(), tr));
            System.out.println(i + ": " + (System.currentTimeMillis() - start));
        }

        System.out.println("Warm JVM: ");
        start = System.currentTimeMillis();
        results.add(Transformation.Util.applyUntilUnchanged(new ExpandTransformation(tr).transform(multiplyAndRenameConflictingDummies(t1, t2)), tr));
        System.out.println("Normal: " + (System.currentTimeMillis() - start));
        for (int i = 1; i < 5; ++i) {
            start = System.currentTimeMillis();
            ParallelExpand pe = new ParallelExpand((Sum) t1, (Sum) t2, i, tr, null);
            results.add(Transformation.Util.applyUntilUnchanged(pe.expand(), tr));
            System.out.println(i + ": " + (System.currentTimeMillis() - start));
        }


        for (int i = 1; i < results.size(); ++i) {
            if (!TensorUtils.equals(results.get(0), results.get(i))) {
                System.out.println(i);
//                System.out.println(results.get(0));
//                System.out.println(results.get(i));
                System.out.println(results.get(0).size());
                System.out.println(results.get(i).size());
                System.out.println("\n\n");
                for (int k = 0; k < results.get(0).size(); ++k) {
                    System.out.println(results.get(0).get(k));
                    if (results.get(i).size() > k)
                        System.out.println(results.get(i).get(k));
                    else System.out.println("-------------------");
                }
            }
            Assert.assertTrue(TensorUtils.equals(results.get(0), results.get(i)));
        }
    }

    @Test
    public void tes2() throws Exception {

        RandomGenerator rw = new Well1024a(new SecureRandom().nextInt());
        for (int i = 0; i < 100; ++i) {
            int seed = rw.nextInt();
            System.out.println(i + "  " + seed);
            CC.resetTensorNames(seed);

            setAntiSymmetric("e_abcd");
//            setAntiSymmetric("D_ac");
//            setAntiSymmetric("B_abc");
//            setAntiSymmetric("A_abc");

            RandomTensor rnd = new RandomTensor(false);
            rnd.reset(seed);

            rnd.addToNamespace(parse("F_a"));
            rnd.addToNamespace(parse("A_ab"));
            rnd.addToNamespace(parse("B_abc"));
            rnd.addToNamespace(parse("D_ac"));
            rnd.addToNamespace(parse("g_ac"));
            rnd.addToNamespace(parse("e_abcd"));


            Sum t1 = (Sum) rnd.nextSum(100, 8, IndicesFactory.EMPTY_INDICES);
            Sum t2 = (Sum) rnd.nextSum(100, 8, IndicesFactory.EMPTY_INDICES);
            t1 = (Sum) parseExpression("e_abcd*e_pqrs = 0").transform(t1);
            t2 = (Sum) parseExpression("e_abcd*e_pqrs = 0").transform(t2);
            System.out.println(t1.size() + "  " + t2.size());
            Transformation tr = new TransformationCollection(
                    EliminateMetricsTransformation.ELIMINATE_METRICS,
                    Tensors.parseExpression("A_ab*B^bac = T^c"),
                    Tensors.parseExpression("A_ab*A^ba = xx"),
                    Tensors.parseExpression("D_ab*D^ba = yy"),
                    EliminateDueSymmetriesTransformation.ELIMINATE_DUE_SYMMETRIES,
                    new LeviCivitaSimplifyTransformation(parseSimple("e_abcd"), true),
                    ExpandAndEliminateTransformation.EXPAND_AND_ELIMINATE,
                    Tensors.parseExpression("A_ab*B^bac = T^c"),
                    Tensors.parseExpression("A_ab*A^ba = xx"),
                    Tensors.parseExpression("D_ab*D^ba = yy"),
                    Tensors.parseExpression("d^a_a = 4")
            );
            t2 = (Sum) ApplyIndexMapping.renameDummy(t2, TensorUtils.getAllIndicesNamesT(t1).toArray());

            Tensor r1 = new ExpandTransformation(tr).transform(multiplyAndRenameConflictingDummies(t1, t2));
            Tensor r3 = new ParallelExpand(t1, t2, 4, tr, null).expand();

            Assert.assertTrue(TensorUtils.isZero(tr.transform(subtract(r1, r3))));
        }
    }


    @Test
    public void tes3() throws Exception {

        long seed = -1721638804;
        CC.resetTensorNames(seed);

//        setAntiSymmetric("e_abcd");
//        setAntiSymmetric("D_ac");
//        setAntiSymmetric("B_abc");
//        setAntiSymmetric("A_abc");

        RandomTensor rnd = new RandomTensor(false);
        rnd.reset(seed);

        rnd.addToNamespace(parse("F_a"));
        rnd.addToNamespace(parse("A_ab"));
        rnd.addToNamespace(parse("B_abc"));
        rnd.addToNamespace(parse("D_ac"));
        rnd.addToNamespace(parse("g_ac"));
        rnd.addToNamespace(parse("e_abcd"));


        Sum t1 = (Sum) rnd.nextSum(80, 6, IndicesFactory.EMPTY_INDICES);
        Sum t2 = (Sum) rnd.nextSum(80, 6, IndicesFactory.EMPTY_INDICES);
        Transformation tr = new TransformationCollection(
                EliminateMetricsTransformation.ELIMINATE_METRICS,
                Tensors.parseExpression("A_ab*B^bac = T^c"),
                Tensors.parseExpression("A_ab*A^ba = xx"),
                Tensors.parseExpression("D_ab*D^ba = yy"),
                EliminateDueSymmetriesTransformation.ELIMINATE_DUE_SYMMETRIES,
                new LeviCivitaSimplifyTransformation(parseSimple("e_abcd"), true),
                ExpandAndEliminateTransformation.EXPAND_AND_ELIMINATE,
                Tensors.parseExpression("A_ab*B^bac = T^c"),
                Tensors.parseExpression("A_ab*A^ba = xx"),
                Tensors.parseExpression("D_ab*D^ba = yy"),
                Tensors.parseExpression("d^a_a = 4")
        );
        t2 = (Sum) ApplyIndexMapping.renameDummy(t2, TensorUtils.getAllIndicesNamesT(t1).toArray());

        Tensor r1 = new ExpandTransformation(tr).transform(multiply(t1, t2));
        Tensor r3 = new ParallelExpand(t1, t2, 1, tr, null).expand();
        Tensor r2 = ExpandUtils.expandPairOfSums(t1, t2, new Transformation[]{tr});

        System.out.println("r12: " + tr.transform(subtract(r1, r2)));
        System.out.println("r13: " + tr.transform(subtract(r1, r3)));
        System.out.println("r23: " + tr.transform(subtract(r2, r3)));

        Assert.assertTrue(TensorUtils.equals(r1, r2));
        Assert.assertTrue(TensorUtils.equals(r1, r3));
        Assert.assertTrue(TensorUtils.equals(r2, r3));
    }

    @Test
    public void test123() throws Exception {
        setAntiSymmetric("e_abcd");
        setAntiSymmetric("D_ac");
        setAntiSymmetric("B_abc");
        setAntiSymmetric("A_abc");

        Tensor t = parse("-16*D_{b}^{i}*B^{b}_{ei}*B_{o}^{qr}*B^{ow}_{q}*B^{e}_{p}^{v}*A_{vw}*A_{r}^{p}-16*D_{b}^{i}*B^{b}_{ei}*B_{o}^{vs}*B^{e}_{ps}*B^{owr}*A_{vw}*A_{r}^{p}+16*D_{b}^{i}*B^{b}_{ei}*B_{t}^{ev}*B_{p}^{tq}*B^{rw}_{q}*A_{r}^{p}*A_{vw}+8*T^{c}*D_{b}^{i}*B^{q}_{co}*B^{bv}_{i}*B^{ow}_{q}*A_{vw}-12*D_{qa}*D_{gr}*F_{t}*F^{b}*e^{cq}_{s}^{r}*A^{a}_{b}*A^{ts}*A_{c}^{g}+24*D_{hc}*B_{q}^{c}_{o}*F^{o}*F_{r}*F_{s}*e_{g}^{q}_{t}^{r}*A^{gh}*A^{ts}+1680*D_{qr}*B^{hf}_{o}*B^{g}_{wh}*F^{o}*F_{t}*e^{wq}_{s}^{r}*A_{gf}*A^{ts}-2*D_{b}^{i}*D^{e}_{s}*B^{b}_{ei}*B_{h}^{d}_{u}*e^{h}_{d}^{r}_{q}*A^{us}*A^{q}_{r}*A^{p}_{p}+12*D_{wr}*F^{c}*F_{c}*F_{t}*F^{d}*e^{wb}_{s}^{r}*A^{ts}*A_{bd}+228*B_{rwe}*F^{e}*F_{q}*F_{t}*F^{c}*F_{b}*e^{wq}_{s}^{r}*A_{c}^{b}*A^{ts}+4*D_{b}^{i}*D^{a}_{s}*B^{b}_{ei}*B_{u}^{rc}*e_{a}^{e}_{cq}*A^{us}*A^{q}_{r}*A^{p}_{p}-24*D_{hc}*B_{q}^{c}_{o}*F^{o}*F_{r}*F_{t}*e_{g}^{q}_{s}^{r}*A^{ts}*A^{gh}+2*D_{b}^{i}*D_{cs}*B^{b}_{qi}*B_{h}^{dc}*e^{h}_{du}^{r}*A^{q}_{r}*A^{us}*A^{p}_{p}+4*D_{b}^{i}*D^{h}_{s}*B^{b}_{ei}*B_{hu}^{c}*e^{re}_{cq}*A^{us}*A^{q}_{r}*A^{p}_{p}-3360*D_{ro}*B^{g}_{qh}*B^{hf}_{w}*F^{o}*F_{t}*e^{wq}_{s}^{r}*A_{gf}*A^{ts}+4*D_{b}^{i}*D_{ds}*B^{b}_{ei}*B^{rdc}*e_{u}^{e}_{cq}*A^{us}*A^{q}_{r}*A^{p}_{p}-156*D^{ac}*D_{oc}*F_{a}*F^{o}*F_{r}*F_{s}*e^{dg}_{t}^{r}*A^{ts}*A_{dg}-8*T^{v}*D_{b}^{i}*B^{b}_{ei}*B^{qe}_{o}*B^{ow}_{q}*A_{vw}+228*B_{qre}*F^{e}*F_{o}*F^{o}*F^{c}*F_{s}*e_{b}^{q}_{t}^{r}*A_{c}^{b}*A^{ts}-12*D_{w}^{a}*B_{r}^{b}_{q}*F_{a}*F_{o}*F^{o}*F_{s}*F_{g}*e^{wq}_{t}^{r}*A_{b}^{g}*A^{ts}+16*D_{b}^{i}*B^{b}_{ei}*B^{ow}_{q}*B^{e}_{p}^{q}*B_{o}^{vr}*A_{vw}*A_{r}^{p}-24*D_{wo}*D_{q}^{e}*F^{o}*F_{t}*e^{wq}_{s}^{d}*A^{ts}*A_{de}*A_{h}^{h}+12*D_{qa}*D_{gw}*F^{c}*F_{s}*e^{wq}_{t}^{b}*A^{a}_{b}*A_{c}^{g}*A^{ts}+2*D_{b}^{i}*D^{g}_{s}*B^{b}_{ui}*B_{h}^{dr}*e^{h}_{dqg}*A^{us}*A^{q}_{r}*A^{p}_{p}-12*D_{qa}*D_{gw}*F_{t}*F^{c}*e^{wq}_{s}^{b}*A^{a}_{b}*A_{c}^{g}*A^{ts}-8*D_{b}^{i}*B^{b}_{ei}*B_{t}^{u}_{s}*B^{ts}_{u}*B^{rwe}*A_{r}^{v}*A_{vw}+8*D_{b}^{i}*B^{b}_{ei}*B^{eu}_{s}*B^{q}_{u}^{s}*B^{rw}_{q}*A_{r}^{v}*A_{vw}-16*D_{b}^{i}*B^{b}_{ei}*B_{t}^{ev}*B^{rtq}*B_{p}^{w}_{q}*A_{r}^{p}*A_{vw}-8*T^{w}*D_{b}^{i}*B^{b}_{ei}*B^{eu}_{s}*B^{v}_{u}^{s}*A_{vw}-24*D_{b}^{i}*B^{b}_{ei}*B^{st}_{o}*B_{t}^{r}_{s}*B^{owe}*A_{r}^{v}*A_{vw}+2*D_{b}^{i}*D^{g}_{s}*B^{b}_{qi}*B_{h}^{d}_{u}*e^{h}_{d}^{r}_{g}*A^{us}*A^{q}_{r}*A^{p}_{p}+2*D_{b}^{i}*D_{as}*B^{b}_{ui}*B_{h}^{d}_{q}*e^{h}_{d}^{ar}*A^{us}*A^{q}_{r}*A^{p}_{p}-12*D_{wr}*D_{q}^{e}*F^{d}*F_{s}*e^{wq}_{t}^{r}*A_{de}*A^{ts}*A_{h}^{h}+8*D_{b}^{i}*B^{b}_{ei}*B_{o}^{qs}*B^{er}_{s}*B^{ow}_{q}*A_{r}^{v}*A_{vw}+12*D_{wr}*D_{q}^{e}*F_{t}*F^{d}*e^{wq}_{s}^{r}*A_{de}*A^{ts}*A_{h}^{h}-2*D_{b}^{i}*D^{g}_{s}*B^{b}_{qi}*B_{h}^{dr}*e^{h}_{dug}*A^{q}_{r}*A^{us}*A^{p}_{p}-12*D_{o}^{e}*D_{rw}*F^{o}*F_{s}*e^{wd}_{t}^{r}*A_{de}*A^{ts}*A_{h}^{h}-16*T_{u}*D_{b}^{i}*B^{vu}_{o}*B^{b}_{ei}*B^{owe}*A_{vw}+24*D_{r}^{a}*B_{w}^{b}_{o}*F_{a}*F^{o}*F_{q}*F_{t}*F_{g}*e^{wq}_{s}^{r}*A_{b}^{g}*A^{ts}+2*D_{b}^{i}*D^{e}_{s}*B^{b}_{ei}*B_{h}^{d}_{q}*e^{h}_{d}^{r}_{u}*A^{q}_{r}*A^{us}*A^{p}_{p}+12*D_{oa}*D_{gr}*F^{o}*F_{s}*e^{bc}_{t}^{r}*A^{a}_{b}*A_{c}^{g}*A^{ts}+24*D_{b}^{i}*B^{bv}_{i}*B^{ts}_{o}*B_{tps}*B^{owr}*A_{vw}*A_{r}^{p}+16*D_{b}^{i}*B^{b}_{ei}*B_{o}^{qs}*B^{ow}_{q}*B^{e}_{ps}*A_{vw}*A^{vp}+16*T^{q}*D_{b}^{i}*B^{ow}_{q}*B^{b}_{ei}*B^{ve}_{o}*A_{vw}-12*D_{go}*D_{ra}*F^{o}*F_{s}*e^{bc}_{t}^{r}*A^{a}_{b}*A_{c}^{g}*A^{ts}+156*D^{ac}*D_{oc}*F_{a}*F^{o}*F_{r}*F_{t}*e^{dg}_{s}^{r}*A^{ts}*A_{dg}-24*D_{wo}*F^{o}*F_{r}*F_{t}*F^{d}*e^{wb}_{s}^{r}*A^{ts}*A_{bd}+4*D_{b}^{i}*D^{h}_{s}*B^{b}_{ei}*B_{hq}^{c}*e_{u}^{e}_{c}^{r}*A^{q}_{r}*A^{us}*A^{p}_{p}-24*D_{r}^{a}*B_{w}^{b}_{o}*F_{a}*F^{o}*F_{q}*F_{s}*F_{g}*e^{wq}_{t}^{r}*A_{b}^{g}*A^{ts}+12*D_{hc}*B_{q}^{c}_{w}*F^{o}*F_{o}*F_{s}*e^{wq}_{tg}*A^{gh}*A^{ts}-1680*D_{qw}*B^{g}_{oh}*B^{hf}_{r}*F^{o}*F_{s}*e^{wq}_{t}^{r}*A_{gf}*A^{ts}-12*D_{hc}*B_{q}^{c}_{w}*F^{o}*F_{o}*F_{t}*e^{wq}_{sg}*A^{ts}*A^{gh}+12*D_{hc}*B_{r}^{c}_{q}*F_{w}*F_{g}*F_{t}*e^{wq}_{s}^{r}*A^{gh}*A^{ts}-2*D_{b}^{i}*D^{g}_{s}*B^{br}_{i}*B_{h}^{d}_{u}*e^{h}_{dqg}*A^{us}*A^{q}_{r}*A^{p}_{p}-12*D_{oa}*D_{gr}*F^{o}*F_{t}*e^{bc}_{s}^{r}*A^{a}_{b}*A^{ts}*A_{c}^{g}+12*D_{go}*D_{ra}*F^{o}*F_{t}*e^{bc}_{s}^{r}*A^{a}_{b}*A^{ts}*A_{c}^{g}-32*D_{b}^{i}*B^{bv}_{i}*B_{t}^{rq}*B_{p}^{t}_{o}*B^{ow}_{q}*A_{r}^{p}*A_{vw}+1680*D_{qw}*B^{g}_{oh}*B^{hf}_{r}*F^{o}*F_{t}*e^{wq}_{s}^{r}*A_{gf}*A^{ts}+12*D_{w}^{a}*B_{r}^{b}_{q}*F_{a}*F_{o}*F^{o}*F_{t}*F_{g}*e^{wq}_{s}^{r}*A_{b}^{g}*A^{ts}-12*D_{rw}*F_{q}*F^{b}*F_{s}*F^{d}*e^{wq}_{t}^{r}*A_{bd}*A^{ts}+2*D_{b}^{i}*D_{cs}*B^{b}_{ui}*B_{h}^{dc}*e^{h}_{d}^{r}_{q}*A^{us}*A^{q}_{r}*A^{p}_{p}+2136*D_{ro}*B_{wqe}*F^{o}*F_{b}*F_{s}*e^{wq}_{t}^{r}*A^{be}*A^{ts}+4*T^{c}*D_{b}^{i}*D^{a}_{s}*B^{b}_{ei}*e_{a}^{e}_{cu}*A^{us}*A^{p}_{p}+24*D_{wo}*D_{q}^{e}*F^{o}*F_{s}*e^{wq}_{t}^{d}*A_{de}*A^{ts}*A_{h}^{h}+2*D_{b}^{i}*D^{g}_{s}*B^{br}_{i}*B_{h}^{d}_{q}*e^{h}_{dug}*A^{q}_{r}*A^{us}*A^{p}_{p}+24*D_{b}^{i}*B^{bv}_{i}*B^{st}_{o}*B_{t}^{r}_{s}*B^{ow}_{p}*A_{vw}*A_{r}^{p}-1704*B_{rwc}*B_{q}^{e}_{o}*F^{o}*F^{c}*F^{h}*F_{s}*e^{wq}_{t}^{r}*A_{he}*A^{ts}-156*D^{ac}*D_{rc}*F_{a}*F^{o}*F_{o}*F_{t}*e^{dg}_{s}^{r}*A^{ts}*A_{dg}-228*B_{rwe}*F^{e}*F_{q}*F^{c}*F_{s}*F_{b}*e^{wq}_{t}^{r}*A_{c}^{b}*A^{ts}-2136*D_{qw}*B_{roe}*F^{o}*F_{b}*F_{s}*e^{wq}_{t}^{r}*A^{be}*A^{ts}-8*T^{w}*D_{b}^{i}*B^{bv}_{i}*B_{t}^{u}_{s}*B^{st}_{u}*A_{vw}-12*D_{hc}*B_{r}^{c}_{q}*F_{w}*F_{g}*F_{s}*e^{wq}_{t}^{r}*A^{gh}*A^{ts}-2136*D_{ro}*B_{wqe}*F^{o}*F_{t}*F_{b}*e^{wq}_{s}^{r}*A^{be}*A^{ts}+4*D_{b}^{i}*D^{a}_{s}*B^{b}_{ei}*B_{qu}^{c}*e_{a}^{e}_{c}^{r}*A^{us}*A^{q}_{r}*A^{p}_{p}+12*D_{qa}*D_{gr}*F_{s}*F^{b}*e^{cq}_{t}^{r}*A^{a}_{b}*A_{c}^{g}*A^{ts}+1704*B_{rwc}*B_{q}^{e}_{o}*F^{o}*F^{c}*F_{t}*F^{h}*e^{wq}_{s}^{r}*A_{he}*A^{ts}+12*D_{rw}*F_{q}*F_{t}*F^{b}*F^{d}*e^{wq}_{s}^{r}*A_{bd}*A^{ts}+3360*D_{ro}*B^{g}_{qh}*B^{hf}_{w}*F^{o}*F_{s}*e^{wq}_{t}^{r}*A_{gf}*A^{ts}-2*D_{b}^{i}*D_{cs}*B^{br}_{i}*B_{h}^{dc}*e^{h}_{duq}*A^{us}*A^{q}_{r}*A^{p}_{p}+12*D_{o}^{e}*D_{rw}*F^{o}*F_{t}*e^{wd}_{s}^{r}*A^{ts}*A_{de}*A_{h}^{h}+156*D^{ac}*D_{rc}*F_{a}*F^{o}*F_{o}*F_{s}*e^{dg}_{t}^{r}*A^{ts}*A_{dg}-1680*D_{qr}*B^{hf}_{o}*B^{g}_{wh}*F^{o}*F_{s}*e^{wq}_{t}^{r}*A_{gf}*A^{ts}+2*D_{b}^{i}*D^{e}_{s}*B^{b}_{ei}*B_{h}^{dr}*e^{h}_{duq}*A^{us}*A^{q}_{r}*A^{p}_{p}+2136*D_{qw}*B_{roe}*F^{o}*F_{b}*F_{t}*e^{wq}_{s}^{r}*A^{be}*A^{ts}");
        t = parseExpression("A_a^a = 0").transform(t);
        t = parseExpression("D_a^a = 0").transform(t);
        t = parseExpression("B_a^ab = 0").transform(t);
        t = parseExpression("e_abcd = 0").transform(t);
        System.out.println(t);
        for (Tensor tensor : t) {
            System.out.println(EliminateDueSymmetriesTransformation.ELIMINATE_DUE_SYMMETRIES.transform(tensor));
        }
    }
}
