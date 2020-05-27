using FluentAssertions;
using System;
using Xunit;

namespace WebApplication1.Tests
{
    public class guid_provider_should
    {
        [Fact]
        public void never_return_a_empty_guid()
        {
            var provider = new GuidProvider();
            var id = provider.Id;

            id.Should().NotBe(Guid.Empty, "Empty guid can't be returned");
        }
    }
}
