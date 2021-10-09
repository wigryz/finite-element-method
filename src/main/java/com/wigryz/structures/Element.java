package com.wigryz.structures;

import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.Setter;
import lombok.ToString;

import java.util.List;

@Getter
@Setter
@AllArgsConstructor
@ToString
public class Element {

    private int id;
    private List<Integer> idList;

    public Element(int id, int id1, int id2, int id3, int id4) {
        this.id = id;
        idList = List.of(id1, id2, id3, id4);
    }
}
