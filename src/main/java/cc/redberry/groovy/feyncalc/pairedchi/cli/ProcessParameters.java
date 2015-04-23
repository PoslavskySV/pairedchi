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
package cc.redberry.groovy.feyncalc.pairedchi.cli;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;

import java.io.File;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public abstract class ProcessParameters {
    @Parameter(names = {"-f", "--force"}, description = "Force overwrite of output files.")
    public Boolean force;

    @Parameter(names = {"-h", "--help"}, description = "output_file")
    public Boolean help;

    public boolean help() {
        return help != null && help;
    }

    public boolean isForceOverwrite() {
        return force != null && force;
    }

    public abstract String getOutputFileName();

    public void validate() {
        if (!isForceOverwrite()) {
            String fileName = getOutputFileName();
            File file = new File(fileName);
            if (file.exists())
                throw new ParameterException("File " + fileName + " already exists. Use -f option to overwrite it.");
        }
    }

    public File getOutputFile() {
        File file = new File(getOutputFileName());
        if (file.exists() && isForceOverwrite()) {
            file.delete();
            file = new File(getOutputFileName());
        }
        return file;
    }
}
