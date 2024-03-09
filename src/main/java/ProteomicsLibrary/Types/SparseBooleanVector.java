/*
 * Copyright 2018-2019 The Hong Kong University of Science and Technology
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

package ProteomicsLibrary.Types;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SparseBooleanVector {

    private Set<Integer> sparseVector;

    public SparseBooleanVector(Set<Integer> sparseVector) {
        this.sparseVector = new HashSet<>(sparseVector);
    }

    public double norm2square() {
        return sparseVector.size();
    }

    public double dot(SparseVector other) {
        double output = 0;
        for (int i : sparseVector) {
            output += other.get(i);
        }
        return output;
    }


    public int dot(SparseBooleanVector other) {
        Set<Integer> intersectedKeys = new HashSet<>(sparseVector);
        intersectedKeys.retainAll(other.sparseVector);
        return intersectedKeys.size();
    }
    public int union(SparseBooleanVector other) {
        Set<Integer> unionKeys = new HashSet<>(sparseVector);
        unionKeys.addAll(other.sparseVector);
        return unionKeys.size();
    }


    public String toString() {
        StringBuilder sb = new StringBuilder(sparseVector.size() * 6);
        for (int idx : sparseVector) {
            sb.append(idx);
            sb.append(";");
        }
        return sb.toString();
    }
}
