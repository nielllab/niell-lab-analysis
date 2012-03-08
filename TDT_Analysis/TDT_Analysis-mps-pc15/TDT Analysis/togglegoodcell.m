function togglegoodcell(src,eventdata)

  global goodcells;

  cellnum = get(src,'UserData');
  goodcells(cellnum) = ~goodcells(cellnum);

  set(src,'FaceAlpha',0.25 - get(src,'FaceAlpha'));
  set(src,'EdgeAlpha',0.25 - get(src,'EdgeAlpha'));
  cellnum