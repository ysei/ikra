class Integer
  def to_bool
    self == 1 ? true : false
  end
end

class TrueClass
  def to_i
    1
  end
end

class FalseClass
  def to_i
    0
  end
end
