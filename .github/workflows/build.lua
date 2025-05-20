local domfilter = require("make4ht-domfilter")

local process = domfilter {
  function(dom)
    for _,mstyle in ipairs(dom:query_selector("mstyle")) do
      local children = mstyle._children
      -- if there is only one child and it is a text node
      -- then we need to change the element to mtext
      -- this will prevent MathJax error
      if #children == 1 and children[1]:is_text() then
        mstyle._name = "mtext"
      end
    end
    return dom
  end

}

Make:match("html", process)
-- for faster compilation, use $ make4ht -m draft -c config.cfg 00-main.tex
if mode == "draft" then
  Make:htlatex {}
else
  Make:htlatex {}
  Make:bibtex {}
  Make:autohtlatex {}
end
