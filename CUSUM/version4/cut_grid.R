grid.point <- within(grid.point, {
  grp.x = cut(x, (0:10)/10, labels = FALSE)
  grp.y = cut(y, (0:10)/10, labels = FALSE)
})
plot(y ~ x, data = grid.point, pch = (15:25)[grp.x], col = grp.y)
abline(v = (1:9)/10)
abline(h = (1:9)/10)