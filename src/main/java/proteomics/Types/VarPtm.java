/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Types;


import proteomics.PTM.InferPTM;

public class VarPtm implements Comparable<VarPtm>{

    public final double mass;
    public final char site;
    public final int priority; // 1 = high; 0 = low.
    public final boolean onlyProteinTerminalIfnc;
    public  String name = null;
    public int position = 0;
    public  String classification = null;
    private final String toString;

    private final int hashCode;
    public double getMass(){return mass;}
    public VarPtm(double mass, char site, int position, String name, String classification, int priority) {
        this.mass = mass;
        this.site = site;
        this.priority = priority;
        this.position = position;
        this.classification = classification;
        this.onlyProteinTerminalIfnc = false;
        this.name = name;
        toString = InferPTM.df3.format(mass) + "@" + site+ "@" + position;
        hashCode = (InferPTM.df3.format(mass) + "@" + site).hashCode();//var mod only differ by mass and site
    }

    public String getStr(){
        return InferPTM.df3.format(mass) + "@" + site;
    }
    public int hashCode() {
        return hashCode;
    }

    @Override
    public String toString() {
        return toString;
    }

    public boolean equals(Object other) {
        if (other instanceof VarPtm) {
            VarPtm temp = (VarPtm) other;
            return temp.hashCode == hashCode;
        } else {
            return false;
        }
    }

    @Override
    public int compareTo(VarPtm other) {
        return 0;
    }
}
