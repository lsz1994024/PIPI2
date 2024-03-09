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



import java.util.Locale;
import java.util.TreeMap;

public class PosMassMap extends TreeMap<Integer, Double> {

    public final int peptideLength;

    public PosMassMap(int peptideLength) {
        super();
        this.peptideLength = peptideLength;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(1000);
        for (Integer co : this.keySet()) {
            sb.append(String.format(Locale.US, "%.3f", this.get(co))); // mass with 3 decimal
            sb.append("@");
            sb.append(co.toString());
            sb.append(";");
        }
        return sb.toString();
    }

    public int hashCode() {
        return this.toString().hashCode();
    }

    public boolean equals(Object other) {
        if (other instanceof PosMassMap) {
            PosMassMap temp = (PosMassMap) other;
            return temp.hashCode() == this.hashCode();
        } else {
            return false;
        }
    }

    public PosMassMap clone() {
        super.clone();
        PosMassMap other = new PosMassMap(peptideLength);
        other.clear();
        for (Integer co : this.keySet()) {
            other.put(co, this.get(co));
        }
        return other;
    }
}