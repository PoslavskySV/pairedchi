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

import cc.redberry.core.context.ContextManager;
import cc.redberry.core.tensor.*;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.utils.TensorUtils;
import cc.redberry.pipe.CUtils;
import cc.redberry.pipe.OutputPort;
import cc.redberry.pipe.Processor;
import cc.redberry.pipe.blocks.ParallelProcessor;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
final class ParallelExpand implements OutputPort<Tensor> {
    private final Sum sum1, sum2;
    private final AtomicLong index = new AtomicLong(0);
    private final AtomicInteger percent = new AtomicInteger(0);
    private final long maxIndex;
    private final int threads;
    private final Transformation transformation;
    private final Log log;

    ParallelExpand(Sum sum1, Sum sum2, int threads, Transformation transformation, Log log) {
        this.sum1 = sum1;
        this.sum2 = (Sum) ApplyIndexMapping.renameDummy(sum2, TensorUtils.getAllIndicesNamesT(sum1).toArray());
        this.threads = threads;
        this.transformation = transformation;
        this.maxIndex = sum1.size() * sum2.size();
        this.log = log;
    }

    @Override
    public Tensor take() {
        final long index = this.index.getAndIncrement();
        if (index >= maxIndex)
            return null;

        if (this.log != null) {
            int p = (int) (100.0 * index / maxIndex);
            int pp = percent.get();
            if (p != pp) {
                if (percent.compareAndSet(pp, p))
                    if (p <= 10 || p % 10 == 0)
                        log.log(percent + "%");
            }
        }

        int i1 = (int) (index / sum2.size());
        int i2 = (int) (index % sum2.size());
        return Tensors.multiply(sum1.get(i1), sum2.get(i2));
    }

    public Tensor expand() {
        ParallelProcessor<Tensor, Tensor> pp = new ParallelProcessor<>(this,
                new Processor<Tensor, Tensor>() {
                    @Override
                    public Tensor process(Tensor input) {
                        return transformation.transform(input);
                    }
                }, 1024, this.threads, ContextManager.getExecutorService());

        SumBuilder sb = new SumBuilder((int) maxIndex);
        for (Tensor tensor : CUtils.it(pp))
            sb.put(transformation.transform(tensor));

        return sb.build();
    }

    public static Tensor parallelExpand(Sum sum1, Sum sum2, Transformation transformation, Log log) {
        return new ParallelExpand(sum1, sum2, 2, transformation, log).expand();
    }

    public interface Log {
        void log(Object o);
    }
}
